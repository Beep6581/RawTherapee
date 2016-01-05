/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2016 Adam Reichold <adam.reichold@t-online.de>
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
#include "threadutils.h"

#include <csignal>
#include <iostream>

#ifdef WIN32
#include <windows.h>
#endif

#if STRICT_MUTEX && !NDEBUG

void MyMutex::checkLock ()
{
    if (locked) {
        std::cerr << "MyMutex already locked!" << std::endl;

#ifdef WIN32
        DebugBreak ();
#else
        raise (SIGTRAP);
#endif
    }

    locked = true;
}

void MyMutex::checkUnlock ()
{
    if (!locked) {
        std::cerr << "MyMutex already unlocked!" << std::endl;

#ifdef WIN32
        DebugBreak ();
#else
        raise (SIGTRAP);
#endif
    }

    locked = false;
}

#endif

#if !TRACE_MYRWMUTEX

void MyReaderLock::acquire ()
{
    if (locked) {
        return;
    }

    Glib::Threads::Mutex::Lock lock (mutex.mutex);

    if (mutex.writerCount == 0) {
        // There's no writer operating, we can increment the writer count which will lock writers.
        ++mutex.writerCount;
    } else if (mutex.readerCount == 0) {
        // The writer count is non null, but a reader can be the owner of the writer lock,
        // which will be the case if the reader count is not zero too.
        while (mutex.writerCount != 0) {
            mutex.cond.wait(mutex.mutex);
        }

        // Then, we can increment the writer count.
        ++mutex.writerCount;
    }

    // Finally, we can increment the reader count as well.
    ++mutex.readerCount;

    locked = true;
}

void MyReaderLock::release ()
{
    if (!locked) {
        return;
    }

    Glib::Threads::Mutex::Lock lock (mutex.mutex);

    // decrement the writer number first...
    --mutex.readerCount;

    if (mutex.readerCount == 0) {
        // ...if no more reader, we decrement the writer count as well...
        --mutex.writerCount;

        // ...and signal the next waiting reader/writer that it's free
        mutex.cond.broadcast ();
    }

    locked = false;
}

void MyWriterLock::acquire ()
{
    if (locked) {
        return;
    }

    Glib::Threads::Mutex::Lock lock (mutex.mutex);

    // The writer count is not zero, so we have to wait for it to be zero again...
    while (mutex.writerCount != 0) {
        mutex.cond.wait (mutex.mutex);
    }

    // ...then we can increment the writer count.
    ++mutex.writerCount;

    locked = true;
}

void MyWriterLock::release ()
{
    if (!locked) {
        return;
    }

    Glib::Threads::Mutex::Lock lock (mutex.mutex);

    // Decrement the writer number first...
    if (--mutex.writerCount == 0) {
        // ...and if the writer count is zero again, we can wake up the next writer or reader.
        mutex.cond.broadcast ();
    }

    locked = false;
}

#else

namespace
{

std::ostream& trace (const char* file, int line)
{
    const auto currentThread = Glib::Threads::Thread::self ();

    return std::cout << currentThread << ":" << file << ":" << line << ": ";
}

}

void MyReaderLock::acquire (const char* file, int line)
{
    if (locked) {
        trace (file, line) << "MyReaderLock is already locked." << std::endl;
        return;
    }

    trace (file, line) << "Acquiring MyReaderLock..." << std::endl;

    Glib::Threads::Mutex::Lock lock (mutex.mutex);

    if (mutex.writerCount == 0) {
        // There's no writer operating, we can increment the writer count which will lock writers.
        ++mutex.writerCount;
    } else if (mutex.readerCount == 0) {
        // The writer count is non null, but a reader can be the owner of the writer lock,
        // which will be the case if the reader count is not zero too.
        while (mutex.writerCount != 0) {
            trace (file, line) << "Waiting for current owner of MyWriterLock..." << std::endl
                               << "\tOwner thread: " << mutex.ownerThread << std::endl
                               << "\tLast writer file: " << mutex.lastWriterFile << std::endl
                               << "\tLast writer line: " << mutex.lastWriterLine << std::endl;

            mutex.cond.wait(mutex.mutex);
        }

        // Then, we can increment the writer count.
        ++mutex.writerCount;

        mutex.ownerThread = Glib::Threads::Thread::self ();
        mutex.lastWriterFile = file;
        mutex.lastWriterLine = line;
    }

    // Finally, we can increment the reader count as well.
    ++mutex.readerCount;

    trace (file, line) << "MyReaderLock is now locked, reader count is " << mutex.readerCount << ", writer count is " << mutex.writerCount << "." << std::endl;
    locked = true;
}

void MyReaderLock::release (const char* file, int line)
{
    if (!locked) {
        trace (file, line) << "MyReaderLock is already unlocked." << std::endl;
        return;
    }

    trace (file, line) << "Releasing MyReaderLock..." << std::endl;

    Glib::Threads::Mutex::Lock lock (mutex.mutex);

    // decrement the writer number first...
    --mutex.readerCount;

    if (mutex.readerCount == 0) {
        // ...if no more reader, we decrement the writer count as well...
        --mutex.writerCount;

        // ...and signal the next waiting reader/writer that it's free
        mutex.cond.broadcast ();

        mutex.ownerThread = nullptr;
        mutex.lastWriterFile = "";
        mutex.lastWriterLine = 0;
    }

    trace (file, line) << "MyReaderLock is now unlocked, reader count is " << mutex.readerCount << ", writer count is " << mutex.writerCount << "." << std::endl;
    locked = false;
}

void MyWriterLock::acquire (const char* file, int line)
{
    if (locked) {
        trace (file, line) << "MyWriterLock is already locked." << std::endl;
        return;
    }

    trace (file, line) << "Acquiring MyWriterLock..." << std::endl;

    Glib::Threads::Mutex::Lock lock (mutex.mutex);

    // The writer count is not zero, so we have to wait for it to be zero again...
    while (mutex.writerCount != 0) {
        trace (file, line) << "Waiting for current owner of MyWriterLock..." << std::endl
                           << "\tOwner thread: " << mutex.ownerThread << std::endl
                           << "\tLast writer file: " << mutex.lastWriterFile << std::endl
                           << "\tLast writer line: " << mutex.lastWriterLine << std::endl;

        mutex.cond.wait (mutex.mutex);
    }

    // ...then we can increment the writer count.
    ++mutex.writerCount;

    mutex.ownerThread = Glib::Threads::Thread::self ();
    mutex.lastWriterFile = file;
    mutex.lastWriterLine = line;

    trace (file, line) << "MyWriterLock is now locked, reader count is " << mutex.readerCount << ", writer count is " << mutex.writerCount << "." << std::endl;
    locked = true;
}

void MyWriterLock::release (const char* file, int line)
{
    if (!locked) {
        trace (file, line) << "MyWriterLock is already unlocked." << std::endl;
        return;
    }

    trace (file, line) << "Releasing MyWriterLock..." << std::endl;

    Glib::Threads::Mutex::Lock lock (mutex.mutex);

    // Decrement the writer number first...
    if (--mutex.writerCount == 0) {
        // ...and if the writer count is zero again, we can wake up the next writer or reader.
        mutex.cond.broadcast ();

        mutex.ownerThread = nullptr;
        mutex.lastWriterFile = "";
        mutex.lastWriterLine = 0;
    }

    trace (file, line) << "MyWriterLock is now unlocked, reader count is " << mutex.readerCount << ", writer count is " << mutex.writerCount << "." << std::endl;
    locked = false;
}

#endif
