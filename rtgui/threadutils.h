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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

// Uncomment this if you want to bypass the CMakeList options and force the values, but do not commit!
//#undef TRACE_MYRWMUTEX
//#define TRACE_MYRWMUTEX 1
//#undef STRICT_MUTEX
//#define STRICT_MUTEX 1

#include <mutex>
#include <thread>
#include <condition_variable>
#include "../rtengine/noncopyable.h"

#if STRICT_MUTEX && NDEBUG
using MyMutexBase = std::mutex;
#else
using MyMutexBase = std::recursive_mutex;
#endif

/**
 * @brief Custom implementation to replace Glib::Threads::Mutex.
 *
 * Glib::Threads::Mutex shows different behaviour on Windows (recursive) and Linux (non-recursive).
 * We therefore use a custom implementation that is optionally recursive and instrumented.
 * It will behave like Glib::Threads::RecMutex (STRICT_MUTEX=0) or Glib::Threads::Mutex (STRICT_MUTEX=1).
 * Debug builds with strict mutexes, will emit a message and crash immediately upon recursive locking.
 */
class MyMutex :
    public rtengine::NonCopyable,
    private MyMutexBase
{
public:
    class MyLock;

    void lock ();
    bool trylock ();
    void unlock ();

#if STRICT_MUTEX && !NDEBUG
    MyMutex();

private:
    bool locked;
    void checkLock ();
    void checkUnlock ();
#endif
};

class MyMutex::MyLock :
    public rtengine::NonCopyable
{
public:
    explicit MyLock (MyMutex& mutex);

    ~MyLock ();

    void acquire ();
    bool try_acquire ();
    void release ();

private:
    MyMutex& mutex;
    bool locked;
};

/**
 * @brief Custom implementation to replace Glib::Threads::RWLock
 */
class MyRWMutex :
    public rtengine::NonCopyable
{
public:
    friend class MyReaderLock;
    friend class MyWriterLock;

    MyRWMutex();

private:
#if TRACE_MYRWMUTEX
    std::thread::id ownerThread;
    const char* lastWriterFile;
    int lastWriterLine;
#endif

    std::mutex mutex;
    std::condition_variable cond;

    std::size_t writerCount;
    std::size_t readerCount;
};

/**
 * @brief Custom implementation to replace Glib::Threads::RWLock::ReaderLock
 */
class MyReaderLock :
    public rtengine::NonCopyable
{
public:
    ~MyReaderLock ();

#if !TRACE_MYRWMUTEX
    explicit MyReaderLock (MyRWMutex& mutex);

    void acquire ();
    void release ();
#else
    explicit MyReaderLock (MyRWMutex& mutex, const char* file, int line);

    void acquire (const char* file, int line);
    void release (const char* file, int line);
#endif

private:
    MyRWMutex& mutex;
    bool locked;
};

/**
 * @brief Custom implementation to replace Glib::Threads::RWLock::WriterLock
 */
class MyWriterLock :
    public rtengine::NonCopyable
{
public:
    ~MyWriterLock ();

#if !TRACE_MYRWMUTEX
    explicit MyWriterLock (MyRWMutex& mutex);

    void acquire ();
    void release ();
#else
    MyWriterLock (MyRWMutex& mutex, const char* file, int line);

    void acquire (const char* file, int line);
    void release (const char* file, int line);
#endif

private:
    MyRWMutex& mutex;
    bool locked;
};

inline void MyMutex::lock ()
{
    MyMutexBase::lock ();

#if STRICT_MUTEX && !NDEBUG
    checkLock ();
#endif
}

inline bool MyMutex::trylock ()
{
    if (MyMutexBase::try_lock ()) {
#if STRICT_MUTEX && !NDEBUG
        checkLock ();
#endif

        return true;
    }

    return false;
}

inline void MyMutex::unlock ()
{
#if STRICT_MUTEX && !NDEBUG
    checkUnlock ();
#endif

    MyMutexBase::unlock ();
}

inline MyMutex::MyLock::MyLock (MyMutex& mutex)
    : mutex (mutex)
    , locked (true)
{
    mutex.lock();
}

inline MyMutex::MyLock::~MyLock ()
{
    if (locked) {
        mutex.unlock ();
    }
}

inline void MyMutex::MyLock::acquire ()
{
    mutex.lock ();
    locked = true;
}
inline bool MyMutex::MyLock::try_acquire ()
{
    return locked = mutex.trylock ();
}

inline void MyMutex::MyLock::release ()
{
    mutex.unlock ();
    locked = false;
}

#if !TRACE_MYRWMUTEX

inline MyReaderLock::MyReaderLock (MyRWMutex& mutex)
    : mutex (mutex)
    , locked (false)
{
    acquire ();
}

inline MyWriterLock::MyWriterLock (MyRWMutex& mutex)
    : mutex (mutex)
    , locked (false)
{
    acquire ();
}

inline MyReaderLock::~MyReaderLock ()
{
    if (locked) {
        release ();
    }
}

inline MyWriterLock::~MyWriterLock ()
{
    if (locked) {
        release ();
    }
}

#else

inline MyReaderLock::MyReaderLock (MyRWMutex& mutex, const char* file, int line)
    : mutex (mutex)
    , locked (false)
{
    acquire (file, line);
}

inline MyWriterLock::MyWriterLock (MyRWMutex& mutex, const char* file, int line)
    : mutex (mutex)
    , locked (false)
{
    acquire (file, line);
}

inline MyReaderLock::~MyReaderLock ()
{
    if (locked) {
        release (__FILE__, __LINE__);
    }
}

inline MyWriterLock::~MyWriterLock ()
{
    if (locked) {
        release (__FILE__, __LINE__);
    }
}

#endif

#if TRACE_MYRWMUTEX
#define MYREADERLOCK(ln, e) MyReaderLock ln(e, __FILE__, __LINE__);
#define MYWRITERLOCK(ln, e) MyWriterLock ln(e, __FILE__, __LINE__);
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
