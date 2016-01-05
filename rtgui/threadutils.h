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

// Uncomment this if you want to bypass the CMakeList options and force the values, but do not commit!
//#undef PROTECT_VECTORS
//#define PROTECT_VECTORS 1
//#undef TRACE_MYRWMUTEX
//#define TRACE_MYRWMUTEX 1
//#undef STRICT_MUTEX
//#define STRICT_MUTEX 1

#include <glibmm/threads.h>

#if STRICT_MUTEX && NDEBUG
using MyMutexBase = Glib::Threads::Mutex;
#else
using MyMutexBase = Glib::Threads::RecMutex;
#endif

/**
 * @brief Custom implementation to replace Glib::Threads::Mutex.
 *
 * Glib::Threads::Mutex shows different behaviour on Windows (recursive) and Linux (non-recursive).
 * We therefore use a custom implementation that is optionally recursive and instrumented.
 * It will behave like Glib::Threads::RecMutex (STRICT_MUTEX=0) or Glib::Threads::Mutex (STRICT_MUTEX=1).
 * Debug builds with strict mutexes, will emit a message and crash immediately upon recursive locking.
 */
class MyMutex : private MyMutexBase
{
public:
    class MyLock;

    MyMutex () = default;
    MyMutex (const MyMutex&) = delete;
    MyMutex& operator= (const MyMutex&) = delete;

    void lock ();
    bool trylock ();
    void unlock ();

#if STRICT_MUTEX && !NDEBUG
private:
    bool locked = false;
    void checkLock ();
    void checkUnlock ();
#endif
};

class MyMutex::MyLock
{
public:
    explicit MyLock (MyMutex& mutex);
    MyLock (MyMutex& mutex, Glib::Threads::NotLock);
    MyLock (MyMutex& mutex, Glib::Threads::TryLock);

    ~MyLock ();

    MyLock (const MyLock&) = delete;
    MyLock& operator= (const MyLock&) = delete;

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
class MyRWMutex
{
public:
    MyRWMutex () = default;
    MyRWMutex (const MyRWMutex&) = delete;
    MyRWMutex& operator= (const MyRWMutex&) = delete;

    friend class MyReaderLock;
    friend class MyWriterLock;

private:
    Glib::Threads::Mutex mutex;
    Glib::Threads::Cond cond;

    std::size_t writerCount = 0;
    std::size_t readerCount = 0;

#if TRACE_MYRWMUTEX
    Glib::Threads::Thread* ownerThread = nullptr;
    const char* lastWriterFile = "";
    int lastWriterLine = 0;
#endif
};

/**
 * @brief Custom implementation to replace Glib::Threads::RWLock::ReaderLock
 */
class MyReaderLock
{
public:
    ~MyReaderLock ();

    MyReaderLock (const MyReaderLock&) = delete;
    MyReaderLock& operator= (const MyReaderLock&) = delete;

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
class MyWriterLock
{
public:
    ~MyWriterLock ();

    MyWriterLock (const MyWriterLock&) = delete;
    MyWriterLock& operator= (const MyWriterLock&) = delete;

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
    if (MyMutexBase::trylock ()) {
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

inline MyMutex::MyLock::MyLock (MyMutex& mutex, Glib::Threads::NotLock)
    : mutex (mutex)
    , locked (false)
{
}

inline MyMutex::MyLock::MyLock (MyMutex& mutex, Glib::Threads::TryLock)
    : mutex (mutex)
    , locked (mutex.trylock ())
{
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
