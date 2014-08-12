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
#include "thumbimageupdater.h"
#include <gtkmm.h>
#include "guiutils.h"
#include "threadutils.h"

#ifdef _OPENMP
#include <omp.h>
#endif 

#define DEBUG(format,args...)
//#define DEBUG(format,args...) printf("ThumbImageUpdate::%s: " format "\n", __FUNCTION__, ## args)

class
ThumbImageUpdater::Impl
{
public:

	struct Job
	{
		Job(ThumbBrowserEntryBase* tbe, bool* priority, bool upgrade,
					ThumbImageUpdateListener* listener):
			tbe_(tbe),
			/*pparams_(pparams),
			height_(height), */
			priority_(priority),
			upgrade_(upgrade),
			listener_(listener)
		{}

		Job():
			tbe_(0),
			priority_(NULL),
			upgrade_(false),
			listener_(0)
		{}

		ThumbBrowserEntryBase* tbe_;
		/*rtengine::procparams::ProcParams pparams_;
		int height_;*/
		bool* priority_;
		bool upgrade_;
		ThumbImageUpdateListener* listener_;
	};

	typedef std::list<Job> JobList;

	Impl():
		active_(0),
		inactive_waiting_(false)
	{
		int threadCount=1;
#if !(defined( WIN32 ) && defined( __x86_64__ ))
		// See Issue 2431 for explanation
		#ifdef _OPENMP
			threadCount=omp_get_num_procs();
		#endif
#endif
		
		threadPool_=new Glib::ThreadPool(threadCount,0);
	}

	Glib::ThreadPool* threadPool_;

	// Need to be a Glib::Threads::Mutex because used in a Glib::Threads::Cond object...
	// This is the only exceptions in RT so far, MyMutex is used everywhere else
	#ifdef WIN32
	Glib::Mutex mutex_;
	#else
	Glib::Threads::Mutex mutex_;
	#endif

	JobList jobs_;

	unsigned int active_;

	bool inactive_waiting_;

	#ifdef WIN32
	Glib::Cond inactive_;
	#else
	Glib::Threads::Cond inactive_;
	#endif

	void
	processNextJob()
	{ 
		Job j;

		{
			#ifdef WIN32
			Glib::Mutex::Lock lock(mutex_);
			#else
			Glib::Threads::Mutex::Lock lock(mutex_);
			#endif

			// nothing to do; could be jobs have been removed
			if ( jobs_.empty() )
			{
				DEBUG("processing: nothing to do (%d)",jobs_.empty());
				return;
			}

			JobList::iterator i;

			// see if any priority jobs exist
			for ( i = jobs_.begin(); i != jobs_.end(); ++i)
			{
				if ( *(i->priority_) )
				{
					DEBUG("processing(priority) %s",i->tbe_->thumbnail->getFileName().c_str());
					break;
				}
			}

			// see if any none upgrade jobs exist
			if ( i == jobs_.end() )
			{
				for ( i = jobs_.begin(); i != jobs_.end(); ++i)
				{
					if ( !i->upgrade_ )
					{
						DEBUG("processing(not-upgrade) %s",i->tbe_->thumbnail->getFileName().c_str());
						break;
					}
				}
			}

			// if none, then use first
			if ( i == jobs_.end() )
			{
				i = jobs_.begin();
				DEBUG("processing(first) %s",i->tbe_->thumbnail->getFileName().c_str());
			}

			// copy found job
			j = *i;

			// remove so not run again
			jobs_.erase(i);
			DEBUG("%d job(s) remaining", int(jobs_.size()) );

			++active_;
		}

		// unlock and do processing; will relock on block exit, then call listener
		double scale = 1.0;
		rtengine::IImage8* img = 0;
        Thumbnail* thm=j.tbe_->thumbnail;

		if ( j.upgrade_ )
		{
			if ( thm->isQuick() )
			{
				img = thm->upgradeThumbImage(thm->getProcParams(), j.tbe_->getPreviewHeight(), scale);
			}
		}
		else
		{
			img = thm->processThumbImage(thm->getProcParams(), j.tbe_->getPreviewHeight(), scale);
		}

		if (img)
		{
			DEBUG("pushing image %s",thm->getFileName().c_str());
			j.listener_->updateImage(img, scale, thm->getProcParams().crop);
		}

		{
			#ifdef WIN32
			Glib::Mutex::Lock lock(mutex_);
			#else
			Glib::Threads::Mutex::Lock lock(mutex_);
			#endif


			if ( --active_ == 0 &&
					inactive_waiting_ )
			{
				inactive_waiting_ = false;
				inactive_.signal();
			}
		}
	}
};

ThumbImageUpdater*
ThumbImageUpdater::getInstance(void)
{
	// this will not be deleted...
	static ThumbImageUpdater* instance_ = 0;
	if ( instance_ == 0 )
	{
		instance_ = new ThumbImageUpdater();
	}
	return instance_;
}

ThumbImageUpdater::ThumbImageUpdater():
	impl_(new Impl())
{
}

void 
ThumbImageUpdater::add(ThumbBrowserEntryBase* tbe, bool* priority, bool upgrade, ThumbImageUpdateListener* l) 
{
	// nobody listening?
	if ( l == 0 )
	{
		return;
	}

	#ifdef WIN32
	Glib::Mutex::Lock lock(impl_->mutex_);
	#else
	Glib::Threads::Mutex::Lock lock(impl_->mutex_);
	#endif

    // look up if an older version is in the queue
	Impl::JobList::iterator i(impl_->jobs_.begin());
	for ( ; i != impl_->jobs_.end(); ++i )
	{
		if ( i->tbe_ == tbe &&
				i->listener_ == l &&
				i->upgrade_ == upgrade )
		{
			DEBUG("updating job %s",tbe->shortname.c_str());
			// we have one, update queue entry, will be picked up by thread when processed
			/*i->pparams_ = params;
			i->height_ = height; */
			i->priority_ = priority;
			return;
		}
	}

	// create a new job and append to queue
	DEBUG("queing job %s",tbe->shortname.c_str());
	impl_->jobs_.push_back(Impl::Job(tbe,priority,upgrade,l));

	DEBUG("adding run request %s",tbe->shortname.c_str());
	impl_->threadPool_->push(sigc::mem_fun(*impl_, &ThumbImageUpdater::Impl::processNextJob));
}


void 
ThumbImageUpdater::removeJobs(ThumbImageUpdateListener* listener)
{
	DEBUG("removeJobs(%p)",listener);

	#ifdef WIN32
	Glib::Mutex::Lock lock(impl_->mutex_);
	#else
	Glib::Threads::Mutex::Lock lock(impl_->mutex_);
	#endif

	for( Impl::JobList::iterator i(impl_->jobs_.begin()); i != impl_->jobs_.end(); )
	{
		if (i->listener_ == listener)
		{
			DEBUG("erasing specific job");
			Impl::JobList::iterator e(i++);
			impl_->jobs_.erase(e);
		}
		else
		{
			++i;
		}
	}

	while ( impl_->active_ != 0 )
	{
		// XXX this is nasty... it would be nicer if we weren't called with
		// this lock held
		GThreadUnLock unlock;
		DEBUG("waiting for running jobs1");
		impl_->inactive_waiting_ = true;
		impl_->inactive_.wait(impl_->mutex_);
	}
}

void 
ThumbImageUpdater::removeAllJobs(void) 
{ 
	DEBUG("stop");

	#ifdef WIN32
	Glib::Mutex::Lock lock(impl_->mutex_);
	#else
	Glib::Threads::Mutex::Lock lock(impl_->mutex_);
	#endif

	impl_->jobs_.clear();

	while ( impl_->active_ != 0 )
	{
		// XXX this is nasty... it would be nicer if we weren't called with
		// this lock held
		GThreadUnLock unlock;
		DEBUG("waiting for running jobs2");
		impl_->inactive_waiting_ = true;
		impl_->inactive_.wait(impl_->mutex_);
	}
}

