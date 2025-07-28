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

#include <glib.h>

#include "rtengine/noncopyable.h"

#include "hidpi.h"

namespace rtengine
{
    class IImage8;

namespace procparams
{

    struct CropParams;

}

}

class ThumbBrowserEntryBase;

class ThumbImageUpdateListener
{
public:
    virtual ~ThumbImageUpdateListener() = default;

    struct ImageUpdate {
        rtengine::IImage8* img;  // New thumbnail image
        hidpi::LogicalSize size;  // Desired logical pixel size
        int device_scale;  // logical to device pixel scaling factor
        double scale;  // scale (??)
        // Why is this a reference? Seems dangerous and could diverge from
        // used values?
        const rtengine::procparams::CropParams& crop;  // crop params used (??)

        ImageUpdate(rtengine::IImage8* a_img, hidpi::LogicalSize a_size,
                    int a_device_scale, double a_scale,
                    const rtengine::procparams::CropParams& a_crop)
                : img(a_img), size(a_size), device_scale(a_device_scale),
                  scale(a_scale), crop(a_crop) {}
    };

    /**
     * @brief Called when thumbnail image is update
     * @param update thumbnail image and associated data
     * @note no locks are held when called back
     */
    virtual void updateImage(ImageUpdate&& update) = 0;
};

class ThumbImageUpdater :
    public rtengine::NonCopyable
{
public:
    /**
     * @brief Singleton entry point.
     *
     * @return Pointer to thumbnail image updater.
     */
    static ThumbImageUpdater* getInstance(void);

    /**
     * @brief Add an thumbnail image update request.
     *
     * Code will add the request to the queue and, if needed, start a pool
     * thread to process it.
     *
     * @param t thumbnail
     * @param priority if \c true then run as soon as possible
     * @param l listener waiting on update
     */
    void add(ThumbBrowserEntryBase* tbe, bool* priority, bool upgrade, bool forceUpgrade, ThumbImageUpdateListener* l);

    /**
     * @brief Remove jobs associated with listener \c l.
     *
     * Jobs being processed will be finished. Will not return till all jobs for
     * \c l have been completed.
     *
     * @param listener jobs associated with this will be stopped
     */
    void removeJobs(ThumbImageUpdateListener* listener);

    /**
     * @brief Stop processing and remove all jobs.
     *
     * Will not return till all running jobs have completed.
     */
    void removeAllJobs(void);

private:

    ThumbImageUpdater();
    ~ThumbImageUpdater();

    class Impl;
    Impl* impl_;
};

/**
 * @brief Singleton boiler plate.
 *
 * To use: \c thumbImageUpdater->start() ,
 */
#define thumbImageUpdater ThumbImageUpdater::getInstance()
