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
#ifndef __GUI_UTILS_
#define __GUI_UTILS_

#include <gtkmm.h>
#include <glibmm.h>
#include "../rtengine/rtengine.h"
#include <sstream>
#include <iostream>

Glib::ustring escapeHtmlChars(const Glib::ustring &src);
bool removeIfThere (Gtk::Container* cont, Gtk::Widget* w, bool increference=true);
void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh);
Glib::ustring removeExtension (const Glib::ustring& filename);
Glib::ustring getExtension (const Glib::ustring& filename);
bool confirmOverwrite (Gtk::Window& parent, const std::string& filename);
void writeFailed (Gtk::Window& parent, const std::string& filename);
void drawCrop (Cairo::RefPtr<Cairo::Context> cr, int imx, int imy, int imw, int imh, int startx, int starty, double scale, const rtengine::procparams::CropParams& cparams);

/** 
 * @brief Lock GTK for critical section.
 *
 * Will unlock on destruction. To use:
 *
 *   <code>
 *     {
 *       GThreadLock lock;
 *       // critical code
 *     }
 *   </code>
 */
class GThreadLock
{
public:
	GThreadLock()
	{
		gdk_threads_enter();
	}
	~GThreadLock()
	{
		gdk_threads_leave();
	}
};

/** 
 * @brief Unlock GTK critical section.
 *
 * Will relock on destruction.
 */
class GThreadUnLock
{
public:
	GThreadUnLock()
	{
		gdk_threads_leave();
	}
	~GThreadUnLock()
	{
		gdk_threads_enter();
	}
};

#ifdef NDEBUG
    // We don't trace mutex
    #undef TRACE_MYRWMUTEX
    #define TRACE_MYRWMUTEX 0
#endif


// Uncomment this if you want to bypass the CMakeList options and force the values
// Of course, DO NOT COMMIT! :)

//#undef PROTECT_VECTORS
//#define PROTECT_VECTORS 1
//#undef TRACE_MYRWMUTEX
//#define TRACE_MYRWMUTEX 1


/**
 * @brief Custom RWLock with debugging feature, to replace the buggy Glib::RWLock (can have negative reader_count value!)
 *
 * It may be slower, but thread safe!
 */
class MyRWMutex {
public:
	Glib::Mutex handlerMutex;
	Glib::Cond  access;
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
class MyReaderLock {

	MyRWMutex& rwMutex;
	bool locked;

	#if TRACE_MYRWMUTEX
	static unsigned int readerLockCounter;
	int locknumber;

public:
	inline MyReaderLock(MyRWMutex& mutex, const char* name, const char* file, const int line) : rwMutex(mutex), locked(false), locknumber(0)
	#else
public:
	inline MyReaderLock(MyRWMutex& mutex) : rwMutex(mutex)
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
		}
		else {
			// The writer count is non null, but we can be the owner of the writer lock
			// It will be the case if the reader count is non null too.
			if (!rwMutex.readerCount) {
				// the mutex is in real write mode, we're waiting to see it null
				#if TRACE_MYRWMUTEX
				std::cout << "  waiting..." << std::endl << "Current writer owner: " << rwMutex.lastWriterFile << " : " << rwMutex.lastWriterLine << std::endl;
				#endif
				while (rwMutex.writerCount)
					rwMutex.access.wait(rwMutex.handlerMutex);
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
			}
			else {
				// The writer count is non null, but a reader can be the owner of the writer lock,
				// it will be the case if the reader count is non null too.
				if (!rwMutex.readerCount) {
					// the mutex is in real write mode, we're waiting to see it null
					#if TRACE_MYRWMUTEX
					std::cout << "  waiting..." << std::endl << "Current writer owner: " << rwMutex.lastWriterFile << " : " << rwMutex.lastWriterLine << std::endl;
					#endif
					while (rwMutex.writerCount)
						rwMutex.access.wait(rwMutex.handlerMutex);
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
		else std::cout << thread << "/" << locknumber << " / already locked by this object - R (lock)" << std::endl;
		#endif
	}
	inline ~MyReaderLock() {
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
		else std::cout << thread << "/" << locknumber << " / already unlocked by this object - R" << std::endl;
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
		else std::cout << thread << "/" << locknumber << " / already unlocked - R (release)" << std::endl;
		#endif
	}
};

/**
 * @brief Custom WriterLock with debugging feature, to replace the buggy Glib::RWLock (can have negative reader_count value!)
 *
 */
class MyWriterLock {

	MyRWMutex& rwMutex;
	bool locked;

	#if TRACE_MYRWMUTEX
	static unsigned int writerLockCounter;
	int locknumber;
public:
	inline MyWriterLock(MyRWMutex& mutex, const char* name, const char* file, const int line) : rwMutex(mutex), locked(false), locknumber(0)
	#else
public:
	inline MyWriterLock(MyRWMutex& mutex) : rwMutex(mutex)
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
			while (rwMutex.writerCount)
				rwMutex.access.wait(rwMutex.handlerMutex);
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
				while (rwMutex.writerCount)
					rwMutex.access.wait(rwMutex.handlerMutex);
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
		else std::cout << thread << "/" << locknumber << " / already locked by this object - W (lock)" << std::endl;
		#endif
	}
	inline ~MyWriterLock() {
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
		else std::cout << thread << "/" << locknumber << " / already unlocked by this object - W" << std::endl;
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
		else std::cout << thread << "/" << locknumber << " / already unlocked by this object - W (release)" << std::endl;
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

/**
 * @brief subclass of Gtk::ScrolledWindow in order to handle the scrollwheel
 */
class MyScrolledWindow : public Gtk::ScrolledWindow {

	bool on_scroll_event (GdkEventScroll* event);

public:
	MyScrolledWindow();
};

/**
 * @brief subclass of Gtk::ComboBox in order to handle the scrollwheel
 */
class MyComboBox : public Gtk::ComboBox {

	bool on_scroll_event (GdkEventScroll* event);

public:
	MyComboBox ();
};

/**
 * @brief subclass of Gtk::ComboBoxText in order to handle the scrollwheel
 */
class MyComboBoxText : public Gtk::ComboBoxText {

	bool on_scroll_event (GdkEventScroll* event);

public:
	MyComboBoxText ();
};

/**
 * @brief subclass of Gtk::SpinButton in order to handle the scrollwheel
 */
class MySpinButton : public Gtk::SpinButton {

protected:
	bool on_scroll_event (GdkEventScroll* event);
	bool on_key_press_event (GdkEventKey* event);

public:
	MySpinButton ();
	void updateSize();
};

/**
 * @brief subclass of Gtk::HScale in order to handle the scrollwheel
 */
class MyHScale : public Gtk::HScale {

	bool on_scroll_event (GdkEventScroll* event);
};

/**
 * @brief subclass of Gtk::FileChooserButton in order to handle the scrollwheel
 */
class MyFileChooserButton : public Gtk::FileChooserButton {

protected:
	bool on_scroll_event (GdkEventScroll* event);

public:
	MyFileChooserButton (const Glib::ustring& title, Gtk::FileChooserAction action=Gtk::FILE_CHOOSER_ACTION_OPEN);
};

/**
 * A class which maintains the last folder for a FileChooserDialog or Button by
 * caching it in a a variable (which can be persisted externally).
 * Each time the user selects a file or folder, the provided variable is updated
 * with the associated folder. The path found in the variable is set in the
 * dialog instance at constructions time of this object.
 */
class FileChooserLastFolderPersister: public Glib::Object {
public:

	/**
	 * Installs this persister on the provided GtkFileChooser instance and
	 * applies the current folder found in @p folderVariable for the dialog.
	 *
	 * @param chooser file chooser to maintain
	 * @param folderVariable variable storage to use for this dialog
	 */
	FileChooserLastFolderPersister(Gtk::FileChooser* chooser, Glib::ustring& folderVariable);

	virtual ~FileChooserLastFolderPersister();

private:

	/**
	 * Signal handler for the GtkFileChooser selection action.
	 */
	void selectionChanged();

	Gtk::FileChooser* chooser;
	Glib::ustring& folderVariable;
	sigc::connection selectionChangedConnetion;

};

typedef enum RTOrientation {
	RTO_Left2Right,
	RTO_Bottom2Top,
	RTO_Right2Left,
	RTO_Top2Bottom
} eRTOrientation;

enum TOITypes {
	TOI_TEXT,
	TOI_ICON
};

typedef enum RTNav {
	NAV_NONE,
	NAV_NEXT,
	NAV_PREVIOUS
} eRTNav;

/**
 * @brief Handle the switch between text and image to be displayed in the HBox (to be used in a button/toolpanel)
 */
class TextOrIcon : public Gtk::HBox {

protected:
	Gtk::Image* imgIcon;
	Gtk::Label* label;
	Glib::ustring filename;
	Glib::ustring labelText;
	Glib::ustring tooltipText;

public:
	TextOrIcon (Glib::ustring filename, Glib::ustring labelTx, Glib::ustring tooltipTx, TOITypes type);
	~TextOrIcon ();

	void switchTo(TOITypes type);
};

/**
 * @brief Define a gradient milestone
 */
class GradientMilestone {
public:
	double position;
	double r;
	double g;
	double b;
	double a;

	GradientMilestone(double _p=0., double _r=0., double _g=0., double _b=0., double _a=0.) {
		position = _p; r = _r; g = _g; b = _b; a = _a;
	}
};

/**
 * @brief Handle point coordinates
 */
template <class T>
class Point {
public:
	T x, y;
	Point() {
		x = T(0);
		y = T(0);
	}

	Point(T coordX, T coordY) {
		x = coordX;
		y = coordY;
	}

	void setCoords(T coordX, T coordY) {
		x = coordX;
		y = coordY;
	}
};

/**
 * @brief Handle backbuffers as automatically as possible
 */
class BackBuffer {

protected:
	int x, y, w, h;  // Rectangle where the colored bar has to be drawn
	Cairo::RefPtr<Cairo::Surface> surface;
	bool dirty;  // mean that the Surface has to be (re)allocated

public:
	BackBuffer();

	// set the destination drawing rectangle; return true if the dimensions are different
	bool setDrawRectangle(Glib::RefPtr<Gdk::Window> window, int newX, int newY, int newW, int newH);

	void copySurface(Glib::RefPtr<Gdk::Window> window, GdkRectangle *rectangle=NULL);
	void copySurface(BackBuffer *destBackBuffer, GdkRectangle *rectangle=NULL);
	void copySurface(Cairo::RefPtr<Cairo::Surface> destSurface, GdkRectangle *rectangle=NULL);

	void setDirty(bool isDirty) { dirty = isDirty; if (!dirty && !surface) dirty = true; }
	bool isDirty() { return dirty; }
	// you have to check if the surface is created thanks to surfaceCreated before starting to draw on it
	bool surfaceCreated() { return surface; }
	Cairo::RefPtr<Cairo::Surface> getSurface() { return surface; }
	void deleteSurface() { surface.clear(); dirty=true; }
	// will let you get a Cairo::Context for Cairo drawing operations
	Cairo::RefPtr<Cairo::Context> getContext() { return Cairo::Context::create(surface); }
	int getWidth() { return w; }
	int getHeight() { return h; }
};


#endif
