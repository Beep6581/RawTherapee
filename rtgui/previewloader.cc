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

#ifdef _OPENMP
#include <omp.h>
#endif

#define DEBUG(format,args...)
//#define DEBUG(format,args...) printf("PreviewLoader::%s: " format "\n", __FUNCTION__, ## args)

class PreviewLoader::Impl
{
public:
    struct Job {
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
    /* Issue 2406
        struct OutputJob
        {
            bool complete;
            int dir_id;
            PreviewLoaderListener* listener;
            FileBrowserEntry* fdn;
        };
    */
    struct JobCompare {
        bool operator()(const Job& lhs, const Job& rhs)
        {
            if ( lhs.dir_id_ == rhs.dir_id_ ) {
                return lhs.dir_entry_ < rhs.dir_entry_;
            }

            return lhs.dir_id_ < rhs.dir_id_;
        }
    };

    typedef std::set<Job, JobCompare> JobSet;

    Impl(): nConcurrentThreads(0)
    {
        int threadCount = 2;
#ifdef _OPENMP
        threadCount = omp_get_num_procs();
#endif

        threadPool_ = new Glib::ThreadPool(threadCount, 0);
    }

    Impl(const Impl&) = delete;

    Glib::ThreadPool* threadPool_;
    MyMutex mutex_;
    JobSet jobs_;
    gint nConcurrentThreads;
// Issue 2406   std::vector<OutputJob *> output_;

    void processNextJob()
    {
        Job j;
// Issue 2406       OutputJob *oj;
        {
            MyMutex::MyLock lock(mutex_);

            // nothing to do; could be jobs have been removed
            if ( jobs_.empty() ) {
                DEBUG("processing: nothing to do");
                return;
            }

            // copy and remove front job
            j = *jobs_.begin();
            jobs_.erase(jobs_.begin());
            DEBUG("processing %s", j.dir_entry_.c_str());
            DEBUG("%d job(s) remaining", jobs_.size());
            /* Issue 2406
                        oj = new OutputJob();
                        oj->complete = false;
                        oj->dir_id = j.dir_id_;
                        oj->listener = j.listener_;
                        oj->fdn = 0;
                        output_.push_back(oj);
            */
        }

        g_atomic_int_inc (&nConcurrentThreads);  // to detect when last thread in pool has run out

        // unlock and do processing; will relock on block exit, then call listener
        // if something got
// Issue 2406       FileBrowserEntry* fdn = 0;
        try {
            Thumbnail* tmb = 0;
            {
                if (Glib::file_test(j.dir_entry_, Glib::FILE_TEST_EXISTS)) {
                    tmb = cacheMgr->getEntry(j.dir_entry_);
                }
            }

            if ( tmb ) {
                j.listener_->previewReady(j.dir_id_, new FileBrowserEntry(tmb, j.dir_entry_));
// Issue 2406               fdn = new FileBrowserEntry(tmb,j.dir_entry_);
            }

        } catch (Glib::Error &e) {} catch(...) {}

        /* Issue 2406
                {
                    // the purpose of the output_ vector is to deliver the previewReady() calls in the same
                    // order as we got the jobs from the jobs_ queue.
                    MyMutex::MyLock lock(mutex_);
                    oj->fdn = fdn;
                    oj->complete = true;
                    while (output_.size() > 0 && output_.front()->complete) {
                        oj = output_.front();
                        if (oj->fdn) {
                            oj->listener->previewReady(oj->dir_id,oj->fdn);
                        }
                        output_.erase(output_.begin());
                        delete oj;
                    }
                }
        */
        bool last = g_atomic_int_dec_and_test (&nConcurrentThreads);

        // signal at end
        if (last && jobs_.empty()) {
            j.listener_->previewsFinished(j.dir_id_);
        }
    }
};

PreviewLoader::PreviewLoader():
    impl_(new Impl())
{
}

PreviewLoader* PreviewLoader::getInstance(void)
{
    static PreviewLoader instance_;
    return &instance_;
}

void PreviewLoader::add(int dir_id, const Glib::ustring& dir_entry, PreviewLoaderListener* l)
{
    // somebody listening?
    if ( l != 0 ) {
        {
            MyMutex::MyLock lock(impl_->mutex_);

            // create a new job and append to queue
            DEBUG("saving job %s", dir_entry.c_str());
            impl_->jobs_.insert(Impl::Job(dir_id, dir_entry, l));
        }

        // queue a run request
        DEBUG("adding run request %s", dir_entry.c_str());
        impl_->threadPool_->push(sigc::mem_fun(*impl_, &PreviewLoader::Impl::processNextJob));
    }
}

void PreviewLoader::removeAllJobs(void)
{
    DEBUG("stop %d", impl_->nConcurrentThreads);
    MyMutex::MyLock lock(impl_->mutex_);
    impl_->jobs_.clear();
}


