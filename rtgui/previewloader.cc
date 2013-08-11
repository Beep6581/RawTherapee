/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <set>
#include "previewloader.h"
#include "guiutils.h"
#include "threadutils.h"
#include "../rtengine/safegtk.h"

#ifdef _OPENMP
#include <omp.h>
#endif 

#define DEBUG(format,args...)
//#define DEBUG(format,args...) printf("PreviewLoader::%s: " format "\n", __FUNCTION__, ## args)

class PreviewLoader::Impl
{
public:
	struct Job
	{
		Job(int dir_id, const Glib::ustring& dir_entry, PreviewLoaderListener* listener):
			dir_id_(dir_id),
			dir_entry_(dir_entry),
			listener_(listener)
		{}

		Job():
			dir_id_(0),
			listener_(0)
		{}

		int dir_id_;
		Glib::ustring dir_entry_;
		PreviewLoaderListener* listener_;
	};

	struct JobCompare
	{
		bool operator()(const Job& lhs, const Job& rhs)
		{
			if ( lhs.dir_id_ == rhs.dir_id_ )
			{
				return lhs.dir_entry_ < rhs.dir_entry_;
			}
			return lhs.dir_id_ < rhs.dir_id_;
		}
	};

	typedef std::set<Job,JobCompare> JobSet;

	Impl():nConcurrentThreads(0)
	{
		int threadCount=2;
		#ifdef _OPENMP
		threadCount=int(omp_get_num_procs()*1.5);
		#endif 
		
		threadPool_=new Glib::ThreadPool(threadCount,0);
	}

	Glib::ThreadPool* threadPool_;
	MyMutex mutex_;
	JobSet jobs_;
	gint nConcurrentThreads;

	void processNextJob()
	{ 
		Job j;
		{
			MyMutex::MyLock lock(mutex_);

			// nothing to do; could be jobs have been removed
			if ( jobs_.empty() )
			{
				DEBUG("processing: nothing to do");
				return;
			}

			// copy and remove front job
			j = *jobs_.begin();
			jobs_.erase(jobs_.begin());
			DEBUG("processing %s",j.dir_entry_.c_str());
			DEBUG("%d job(s) remaining",jobs_.size());
		}

		g_atomic_int_inc (&nConcurrentThreads);  // to detect when last thread in pool has run out

		// unlock and do processing; will relock on block exit, then call listener
		// if something got
		try {
			Thumbnail* tmb = 0;
			{
				if (safe_file_test(j.dir_entry_, Glib::FILE_TEST_EXISTS))
				{
					tmb = cacheMgr->getEntry(j.dir_entry_);
				}
			}

			// we got something, so notify listener
			if ( tmb )
			{
				j.listener_->previewReady(j.dir_id_,new FileBrowserEntry(tmb,j.dir_entry_));
			}
		} catch (Glib::Error &e){} catch(...){}

		bool last = g_atomic_int_dec_and_test (&nConcurrentThreads);

		// signal at end
		if (last && jobs_.empty()) j.listener_->previewsFinished(j.dir_id_);
	}
};

PreviewLoader::PreviewLoader():
	impl_(new Impl())
{
}

PreviewLoader* PreviewLoader::getInstance(void)
{
	// this will not be deleted...
	static PreviewLoader* instance_ = NULL;
	if ( instance_ == NULL )
	{
        static MyMutex smutex_;
        MyMutex::MyLock lock(smutex_);

        if ( instance_ == NULL ) instance_ = new PreviewLoader();
	}

	return instance_;
}

void PreviewLoader::add(int dir_id, const Glib::ustring& dir_entry, PreviewLoaderListener* l)
{
	// somebody listening?
	if ( l != 0 )
	{
        {
            MyMutex::MyLock lock(impl_->mutex_);

            // create a new job and append to queue
            DEBUG("saving job %s",dir_entry.c_str());
            impl_->jobs_.insert(Impl::Job(dir_id,dir_entry,l));
        }

		// queue a run request
		DEBUG("adding run request %s",dir_entry.c_str());
		impl_->threadPool_->push(sigc::mem_fun(*impl_, &PreviewLoader::Impl::processNextJob));
	}
}

void PreviewLoader::removeAllJobs(void) 
{ 
	DEBUG("stop %d",impl_->nConcurrentThreads);
	MyMutex::MyLock lock(impl_->mutex_);
	impl_->jobs_.clear();
}


