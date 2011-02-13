/*
 * rtcommon.h
 *
 *  Created on: Jan 19, 2011
 *      Author: gabor
 */

#ifndef RTCOMMON_H_
#define RTCOMMON_H_

#include "config.h"
#ifdef QTBUILD
#include <QString>
#include <QList>
#include <QStringList>
#include <QMutex>
#include <QWaitCondition>
#include <QImage>
namespace rtengine {
typedef QString String;
typedef QVector<float> FloatVector;
typedef QVector<int> IntVector;
typedef QStringList StringList;
typedef QImage DisplayImage;
typedef QMutex Mutex;
typedef QWaitCondition Condition;

#define String2PChar(a)  (a.toAscii().constData())
#define String2StdString(a)  (a.toStdString())
#define StdString2String(a)  (QString(a.c_str()))
}

#else
#include <glibmm.h>
#include <cairomm/cairomm.h>
#include <vector>
namespace rtengine {
typedef Glib::ustring String;
typedef std::vector<float> FloatVector;
typedef std::vector<int> 	IntVector;
typedef std::vector<String> StringList;
typedef Cairo::RefPtr<Cairo::ImageSurface> DisplayImage;
typedef Glib::Mutex	Mutex;
typedef Glib::Cond Condition;

#define String2PChar(a)  (a.c_str())
#define String2StdString(a)  (a)
#define StdString2String(a)  (a)
}
#endif

#endif /* RTCOMMON_H_ */
