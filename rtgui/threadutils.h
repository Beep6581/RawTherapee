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


#ifndef _THREADUTILS_
#define _THREADUTILS_

#include <glibmm.h>
#include <csignal>  // for raise()
#include <iostream>

#ifdef WIN32
#include <windows.h>
#endif


#ifdef NDEBUG
// We don't trace mutex
#undef TRACE_MYRWMUTEX
#define TRACE_MYRWMUTEX 0
#endif


// Uncomment this if you want to bypass the CMakeList options and force the values
// Of course, DO NOT COMMIT!

//#undef PROTECT_VECTORS
//#define PROTECT_VECTORS 1
//#undef TRACE_MYRWMUTEX
//#define TRACE_MYRWMUTEX 1
//#undef STRICT_MUTEX
//#define STRICT_MUTEX 1

/**
 * @brief Custom Mutex to replace Glib::Threads::Mutex, which behave differently on windows (recursive) and linux (non-recursive), by a recursive and "debugable" one
 *
 * This implementation will behave like a Glib::Threads::RecMutex (STRICT_MUTEX=0) or a Glib::Threads::Mutex (STRICT_MUTEX=1), but in this case, the application will
 * crash instead of freezing.
 *
 * In Debug builds, a printf will let you know that the MyMutex was already locked
 *
 * The default and recommended mode is STRICT_MUTEX=1
 */

#ifdef WIN32
class MyMutex : public Glib::RecMutex
{
#else
class MyMutex : public Glib::Threads::RecMutex
{
#endif

#if STRICT_MUTEX || !defined(NDEBUG)
private:
    bool alreadyLocked;
#endif

public:
    class MyLock;

#if STRICT_MUTEX || !defined(NDEBUG)
    MyMutex() : alreadyLocked(false) {}
#else
    MyMutex() {}
#endif

    void lock()
    {
#ifdef WIN32
        Glib::RecMutex::lock();
#else
        Glib::Threads::RecMutex::lock();
#endif
#if STRICT_MUTEX || !defined(NDEBUG)

        if (alreadyLocked) {
#ifndef NDEBUG
            std::cout << "Warning: MyMutex already locked!" << std::endl; // breakpoint
#endif
#if STRICT_MUTEX
#ifndef NDEBUG
#ifdef WIN32
            DebugBreak();
#else
            raise(SIGTRAP);
#endif
#else
            raise(SIGINT);
#endif
#endif
        }

        alreadyLocked = true;
#endif
    }

    bool trylock()
    {
#ifdef WIN32

        if (Glib::RecMutex::trylock())
#else
        if (Glib::Threads::RecMutex::trylock())
#endif
        {
#if STRICT_MUTEX || !defined(NDEBUG)

            if (alreadyLocked) {
#ifndef NDEBUG
                std::cout << "Warning: MyMutex already locked!" << std::endl; // breakpoint
#endif
#if STRICT_MUTEX
#ifndef NDEBUG
#ifdef WIN32
                DebugBreak();
#else
                raise(SIGTRAP);
#endif
#else
                raise(SIGINT);
#endif
#endif
            }

            alreadyLocked = true;
#endif
            return true;
        }

        return false;
    }

    // Warning: the base class of MyMutex is RecMutex, but the mutex is said "unlocked" on first occurrence of "unlock", to avoid overhead.
    void unlock()
    {
#if STRICT_MUTEX || !defined(NDEBUG)
        alreadyLocked = false;
#endif
#ifdef WIN32
        Glib::RecMutex::unlock();
#else
        Glib::Threads::RecMutex::unlock();
#endif
    }
};


// Class copied from the Glibmm source code, to provide a workaround of the behavior's difference between Linux and Windows
class MyMutex::MyLock
{
public:
    explicit inline MyLock(MyMutex& mutex) : mutex_  (mutex), locked_ (true)
    {
        mutex_.lock();
    }
#ifdef WIN32
    inline MyLock(MyMutex& mutex, Glib::NotLock) : mutex_  (mutex), locked_ (false) {}
    inline MyLock(MyMutex& mutex, Glib::TryLock) : mutex_  (mutex), locked_ (mutex.trylock()) {}
#else
    inline MyLock(MyMutex& mutex, Glib::Threads::NotLock) : mutex_  (mutex), locked_ (false) {}
    inline MyLock(MyMutex& mutex, Glib::Threads::TryLock) : mutex_  (mutex), locked_ (mutex.trylock()) {}
#endif
    inline ~MyLock()
    {
        if(locked_) {
            mutex_.unlock();
        }
    }

    inline void acquire()
    {
        mutex_.lock();
        locked_ = true;
    }
    inline bool try_acquire()
    {
        locked_ = mutex_.trylock();
        return locked_;
    }
    inline void release()
    {
        mutex_.unlock();
        locked_ = false;
    }
    inline bool locked() const
    {
        return locked_;
    }

private:
    MyMutex& mutex_;
    bool     locked_;

    // noncopyable
    MyLock(const MyMutex::Lock&);
    MyMutex::Lock& operator=(const MyMutex::Lock&);
};


/**
 * @brief Custom RWLock with debugging feature, to replace the buggy Glib::RWLock (can have negative reader_count value!)
 *
 * It may be slower, but thread safe!
 */
class MyRWMutex
{
public:
#ifdef WIN32
    Glib::Mutex handlerMutex;  // Having a recursive or non-recursive mutex is not important here, so we can use Glib::Mutex
    Glib::Cond  access;
#else
    Glib::Threads::Mutex handlerMutex;  // Having a recursive or non-recursive mutex is not important here, so we can use Glib::Threads::Mutex
    Glib::Threads::Cond  access;
#endif
    size_t      writerCount;
    size_t      readerCount;
#if TRACE_MYRWMUTEX
    Glib::ustring lastWriterFile;
    int           lastWriterLine;
    // Unfortunately, ownerThread may not be the culprit of a deadlock, it can be another concurrent Reader...
    void*         ownerThread;

    MyRWMutex() : writerCount(0), readerCount(0), lastWriterLine(0), ownerThread(NULL) {}
#else
    MyRWMutex() : writerCount(0), readerCount(0) {}
#endif
};

/**
 * @brief Custom ReaderLock with debugging feature, to replace the buggy Glib::RWLock (can have negative reader_count value!)
 *
 */
class MyReaderLock
{

    MyRWMutex& rwMutex;
    bool locked;

#if TRACE_MYRWMUTEX
    static unsigned int readerLockCounter;
    int locknumber;

public:
    inline MyReaderLock(MyRWMutex& mutex, const char* name, const char* file, const int line) : rwMutex(mutex), locked(false), locknumber(0)
#else
public:
    inline MyReaderLock(MyRWMutex & mutex) : rwMutex(mutex)
#endif

    {
        // to operate safely
        rwMutex.handlerMutex.lock();

#if TRACE_MYRWMUTEX
        locknumber = readerLockCounter++;
        void* thread = Glib::Thread::self();
        std::cout << thread << "/" << locknumber << ":" << name <<  " / " << file << " : " << line << " - locking - R";
#endif

        if (!rwMutex.writerCount) {
            // There's no writer operating, we can increment the writer count which will lock writers
            ++rwMutex.writerCount;
#if TRACE_MYRWMUTEX
            std::cout << " ++ new owner";
#endif
        } else {
            // The writer count is non null, but we can be the owner of the writer lock
            // It will be the case if the reader count is non null too.
            if (!rwMutex.readerCount) {
                // the mutex is in real write mode, we're waiting to see it null
#if TRACE_MYRWMUTEX
                std::cout << "  waiting..." << std::endl << "Current writer owner: " << rwMutex.lastWriterFile << " : " << rwMutex.lastWriterLine << std::endl;
#endif

                while (rwMutex.writerCount) {
                    rwMutex.access.wait(rwMutex.handlerMutex);
                }

                ++rwMutex.writerCount;
#if TRACE_MYRWMUTEX
                rwMutex.lastWriterFile = file;
                rwMutex.lastWriterLine = line;
                rwMutex.ownerThread = thread;
                std::cout << thread << "/" << locknumber << ":" << name <<  " / " << file << " : " << line << " - locking - R ++ new owner";
#endif
            }
        }

        // then we can increment the reader count
        ++rwMutex.readerCount;

#if TRACE_MYRWMUTEX
        std::cout << " - ReaderCount: " << rwMutex.readerCount << " - WriterCount: " << rwMutex.writerCount << std::endl;
#endif

        rwMutex.handlerMutex.unlock();

        locked = true;
    }
#if TRACE_MYRWMUTEX
    // locks the MyRWMutex with Read access if this MyReaderLock has not already locked it, otherwise return safely
    inline void acquire(const char* file, const int line)
#else
    // locks the MyRWMutex with Read access if this MyReaderLock has not already locked it, otherwise return safely
    inline void acquire()
#endif
    {
#if TRACE_MYRWMUTEX
        void* thread = Glib::Thread::self();
#endif

        if (!locked) {
            // to operate safely
            rwMutex.handlerMutex.lock();

#if TRACE_MYRWMUTEX
            std::cout << thread << "/" << locknumber << ":" << file << " : " << line << " - locking - R (lock)";
#endif

            if (!rwMutex.writerCount) {
                // There's no writer operating, we can increment the writer count which will lock writers
                ++rwMutex.writerCount;
#if TRACE_MYRWMUTEX
                std::cout << " ++ new owner";
#endif
            } else {
                // The writer count is non null, but a reader can be the owner of the writer lock,
                // it will be the case if the reader count is non null too.
                if (!rwMutex.readerCount) {
                    // the mutex is in real write mode, we're waiting to see it null
#if TRACE_MYRWMUTEX
                    std::cout << "  waiting..." << std::endl << "Current writer owner: " << rwMutex.lastWriterFile << " : " << rwMutex.lastWriterLine << std::endl;
#endif

                    while (rwMutex.writerCount) {
                        rwMutex.access.wait(rwMutex.handlerMutex);
                    }

                    ++rwMutex.writerCount;
#if TRACE_MYRWMUTEX
                    rwMutex.lastWriterFile = file;
                    rwMutex.lastWriterLine = line;
                    rwMutex.ownerThread = thread;
                    std::cout << thread << "/" << locknumber << ":" << file << " : " << line << " - locking - R (lock) ++ new owner";
#endif
                }
            }

            // then we can increment the reader count
            ++rwMutex.readerCount;

#if TRACE_MYRWMUTEX
            std::cout << " - ReaderCount: " << rwMutex.readerCount << " - WriterCount: " << rwMutex.writerCount << std::endl;
#endif

            rwMutex.handlerMutex.unlock();

            locked = true;
        }

#if TRACE_MYRWMUTEX
        else {
            std::cout << thread << "/" << locknumber << " / already locked by this object - R (lock)" << std::endl;
        }

#endif
    }
    inline ~MyReaderLock()
    {
#if TRACE_MYRWMUTEX
        void* thread = Glib::Thread::self();
#endif

        if (locked) {
            // to operate safely
            rwMutex.handlerMutex.lock();

            // decrement the writer number first
            --rwMutex.readerCount;

#if TRACE_MYRWMUTEX
            std::cout << thread << "/" << locknumber << " / unlocking - R - ReaderCount: " << rwMutex.readerCount;
#endif

            if (!rwMutex.readerCount) {
                // no more reader, so we decrement the writer count
                --rwMutex.writerCount;
#if TRACE_MYRWMUTEX
                rwMutex.lastWriterFile = "";
                rwMutex.lastWriterLine = 0;
                rwMutex.ownerThread = NULL;
                std::cout << " -- new owner possible!" << " >>> ReaderCount: " << rwMutex.readerCount << " - WriterCount: " << rwMutex.writerCount;
#endif
                // and signal the next waiting reader/writer that it's free
                rwMutex.access.broadcast();
            }

#if TRACE_MYRWMUTEX
            std::cout << std::endl;
#endif

            rwMutex.handlerMutex.unlock();
        }

#if TRACE_MYRWMUTEX
        else {
            std::cout << thread << "/" << locknumber << " / already unlocked by this object - R" << std::endl;
        }

#endif
    }
#if TRACE_MYRWMUTEX
    // releases the MyRWMutex with Write access if this MyWriterLock has already locked it, otherwise return safely
    inline void release(const char* file, const int line)
#else
    // releases the MyRWMutex with Write access if this MyWriterLock has already locked it, otherwise return safely
    inline void release()
#endif
    {
#if TRACE_MYRWMUTEX
        void* thread = Glib::Thread::self();
#endif

        if (locked) {
            // to operate safely
            rwMutex.handlerMutex.lock();

            // decrement the writer number first
            --rwMutex.readerCount;

#if TRACE_MYRWMUTEX
            std::cout << thread << "/" << locknumber << " / unlocking - R (release) - ReaderCount: " << rwMutex.readerCount;
#endif

            if (!rwMutex.readerCount) {
                // no more reader, so we decrement the writer count
                --rwMutex.writerCount;
#if TRACE_MYRWMUTEX
                rwMutex.lastWriterFile = "";
                rwMutex.lastWriterLine = 0;
                rwMutex.ownerThread = NULL;
                std::cout << " -- new owner possible!" << " >>> ReaderCount: " << rwMutex.readerCount << " - WriterCount: " << rwMutex.writerCount;
#endif
                // and signal the next waiting reader/writer that it's free
                rwMutex.access.broadcast();
            }

#if TRACE_MYRWMUTEX
            std::cout << std::endl;
#endif

            rwMutex.handlerMutex.unlock();

            locked = false;
        }

#if TRACE_MYRWMUTEX
        else {
            std::cout << thread << "/" << locknumber << " / already unlocked - R (release)" << std::endl;
        }

#endif
    }
};

/**
 * @brief Custom WriterLock with debugging feature, to replace the buggy Glib::RWLock (can have negative reader_count value!)
 *
 */
class MyWriterLock
{

    MyRWMutex& rwMutex;
    bool locked;

#if TRACE_MYRWMUTEX
    static unsigned int writerLockCounter;
    int locknumber;
public:
    inline MyWriterLock(MyRWMutex& mutex, const char* name, const char* file, const int line) : rwMutex(mutex), locked(false), locknumber(0)
#else
public:
    inline MyWriterLock(MyRWMutex & mutex) : rwMutex(mutex)
#endif
    {
        // to operate safely
        rwMutex.handlerMutex.lock();

#if TRACE_MYRWMUTEX
        locknumber = writerLockCounter++;
        void* thread = Glib::Thread::self();
        std::cout << thread << "/" << locknumber << ":" << name <<  " / " << file << " : " << line << " - locking - W";
#endif

        if (rwMutex.writerCount) {
            // The writer count is non null, so we have to wait for it to be null again
#if TRACE_MYRWMUTEX
            std::cout << "  waiting..." << std::endl << "Current writer owner: " << rwMutex.lastWriterFile << " : " << rwMutex.lastWriterLine << std::endl;
#endif

            while (rwMutex.writerCount) {
                rwMutex.access.wait(rwMutex.handlerMutex);
            }

#if TRACE_MYRWMUTEX
            std::cout << thread << "/" << locknumber << ":" << file << " : " << line << " - locking - W";
#endif
        }

        // then we can increment the writer count
        ++rwMutex.writerCount;

#if TRACE_MYRWMUTEX
        rwMutex.lastWriterFile = file;
        rwMutex.lastWriterLine = line;
        rwMutex.ownerThread = thread;
        std::cout << " ++ new owner <<< ReaderCount: " << rwMutex.readerCount << " - WriterCount: " << rwMutex.writerCount << std::endl;
#endif

        rwMutex.handlerMutex.unlock();

        locked = true;
    }
#if TRACE_MYRWMUTEX
    // locks the MyRWMutex with Read access if this MyReaderLock has not already locked it, otherwise return safely
    inline void acquire(const char* file, const int line)
#else
    // locks the MyRWMutex with Read access if this MyReaderLock has not already locked it, otherwise return safely
    inline void acquire()
#endif
    {
#if TRACE_MYRWMUTEX
        void* thread = Glib::Thread::self();
#endif

        if (!locked) {
            // to operate safely
            rwMutex.handlerMutex.lock();

#if TRACE_MYRWMUTEX
            std::cout << thread << "/" << locknumber << ":" << file << " : " << line << " - locking - W (lock)";
#endif

            if (rwMutex.writerCount) {
                // The writer count is non null, so we have to wait for it to be null again
#if TRACE_MYRWMUTEX
                std::cout << "  waiting..." << std::endl << "Current writer owner: " << rwMutex.lastWriterFile << " : " << rwMutex.lastWriterLine << std::endl;
#endif

                while (rwMutex.writerCount) {
                    rwMutex.access.wait(rwMutex.handlerMutex);
                }

#if TRACE_MYRWMUTEX
                std::cout << thread << "/" << locknumber << ":" << file << " : " << line << " - locking - W (lock)";
#endif
            }

            // then we can increment the reader count
            ++rwMutex.writerCount;

#if TRACE_MYRWMUTEX
            rwMutex.lastWriterFile = file;
            rwMutex.lastWriterLine = line;
            rwMutex.ownerThread = thread;
            std::cout << " ++ new owner <<< ReaderCount: " << rwMutex.readerCount << " - WriterCount: " << rwMutex.writerCount << std::endl;
#endif

            rwMutex.handlerMutex.unlock();

            locked = true;
        }

#if TRACE_MYRWMUTEX
        else {
            std::cout << thread << "/" << locknumber << " / already locked by this object - W (lock)" << std::endl;
        }

#endif
    }
    inline ~MyWriterLock()
    {
#if TRACE_MYRWMUTEX
        void* thread = Glib::Thread::self();
#endif

        if (locked) {
            // to operate safely
            rwMutex.handlerMutex.lock();

            // decrement the writer number first
            --rwMutex.writerCount;

#if TRACE_MYRWMUTEX
            std::cout << thread << "/" << locknumber << " / unlocking - W";
#endif

            if (!rwMutex.writerCount) {
#if TRACE_MYRWMUTEX
                rwMutex.lastWriterFile = "";
                rwMutex.lastWriterLine = 0;
                rwMutex.ownerThread = NULL;
                std::cout << " -- new owner possible!";
#endif
                // The writer count is null again, so we can wake up the next writer or reader
                rwMutex.access.broadcast();
            }

#if TRACE_MYRWMUTEX
            std::cout << " <<< ReaderCount: " << rwMutex.readerCount << " - WriterCount: " << rwMutex.writerCount << std::endl;
#endif

            rwMutex.handlerMutex.unlock();
        }

#if TRACE_MYRWMUTEX
        else {
            std::cout << thread << "/" << locknumber << " / already unlocked by this object - W" << std::endl;
        }

#endif
    }
#if TRACE_MYRWMUTEX
    // releases the MyRWMutex with Write access if this MyWriterLock has already locked it, otherwise return safely
    inline void release(const char* file, const int line)
#else
    // releases the MyRWMutex with Write access if this MyWriterLock has already locked it, otherwise return safely
    inline void release()
#endif
    {
#if TRACE_MYRWMUTEX
        void* thread = Glib::Thread::self();
#endif

        if (locked) {
            // to operate safely
            rwMutex.handlerMutex.lock();

            // decrement the writer number first
            --rwMutex.writerCount;

#if TRACE_MYRWMUTEX
            std::cout << thread << "/" << locknumber << " / unlocking - W (release)";
#endif

            if (!rwMutex.writerCount) {
#if TRACE_MYRWMUTEX
                rwMutex.lastWriterFile = "";
                rwMutex.lastWriterLine = 0;
                rwMutex.ownerThread = NULL;
                std::cout << " -- new owner possible!";
#endif
                // The writer count is null again, so we can wake up the next writer or reader
                rwMutex.access.broadcast();
            }

#if TRACE_MYRWMUTEX
            std::cout << " <<< ReaderCount: " << rwMutex.readerCount << " - WriterCount: " << rwMutex.writerCount << std::endl;
#endif

            rwMutex.handlerMutex.unlock();

            locked = false;
        }

#if TRACE_MYRWMUTEX
        else {
            std::cout << thread << "/" << locknumber << " / already unlocked by this object - W (release)" << std::endl;
        }

#endif
    }
};

#if TRACE_MYRWMUTEX
#define MYREADERLOCK(ln, e) MyReaderLock ln(e, #e, __FILE__, __LINE__);
#define MYWRITERLOCK(ln, e) MyWriterLock ln(e, #e, __FILE__, __LINE__);
#define MYREADERLOCK_ACQUIRE(ln) ln.acquire(__FILE__, __LINE__);
#define MYWRITERLOCK_ACQUIRE(ln) ln.acquire(__FILE__, __LINE__);
#define MYREADERLOCK_RELEASE(ln) ln.release(__FILE__, __LINE__);
#define MYWRITERLOCK_RELEASE(ln) ln.release(__FILE__, __LINE__);
#else
#define MYREADERLOCK(ln, e) MyReaderLock ln(e);
#define MYWRITERLOCK(ln, e) MyWriterLock ln(e);
#define MYREADERLOCK_ACQUIRE(ln) ln.acquire();
#define MYWRITERLOCK_ACQUIRE(ln) ln.acquire();
#define MYREADERLOCK_RELEASE(ln) ln.release();
#define MYWRITERLOCK_RELEASE(ln) ln.release();
#endif

#ifdef PROTECT_VECTORS
#define IFPV_MYREADERLOCK(l, e) MYREADERLOCK(l, e)
#define IFPV_MYWRITERLOCK(l, e) MYWRITERLOCK(l, e)
#define IFPV_MYREADERLOCK_ACQUIRE(l) MYREADERLOCK_ACQUIRE(l)
#define IFPV_MYWRITERLOCK_ACQUIRE(l) MYWRITERLOCK_ACQUIRE(l)
#define IFPV_MYREADERLOCK_RELEASE(l) MYREADERLOCK_RELEASE(l)
#define IFPV_MYWRITERLOCK_RELEASE(l) MYWRITERLOCK_RELEASE(l)
#else
#define IFPV_MYREADERLOCK(l, e)
#define IFPV_MYWRITERLOCK(l, e)
#define IFPV_MYREADERLOCK_ACQUIRE(l)
#define IFPV_MYWRITERLOCK_ACQUIRE(l)
#define IFPV_MYREADERLOCK_RELEASE(l)
#define IFPV_MYWRITERLOCK_RELEASE(l)
#endif

#endif /* _THREADUTILS_ */
