#include <cstdio>
#include <jpeglib.h>
#include <jerror.h>
#include "jpeg.h"

/*
 * jdatasrc.c
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * Modified 2009-2010 by Guido Vollbeding.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains decompression data source routines for the case of
 * reading JPEG data from a file (or any stdio stream).  While these routines
 * are sufficient for most applications, some will want to use a different
 * source manager.
 * IMPORTANT: we assume that fread() will correctly transcribe an array of
 * JOCTETs from 8-bit-wide elements on external storage.  If char is wider
 * than 8 bits on your machine, you may need to do some tweaking.
 */

/* this is not a core library module, so it doesn't define JPEG_INTERNALS */
//#include "jinclude.h"

#define JFREAD(file,buf,sizeofbuf)  \
  ((size_t) fread((void *) (buf), (size_t) 1, (size_t) (sizeofbuf), (file)))
#define JFWRITE(file,buf,sizeofbuf)  \
  ((size_t) fwrite((const void *) (buf), (size_t) 1, (size_t) (sizeofbuf), (file)))




/* Expanded data source object for stdio input */
namespace
{

typedef struct {
    struct jpeg_source_mgr pub;   /* public fields */
    jmp_buf error_jmp_buf; /* error handler for this instance */

    FILE * infile;        /* source stream */
    JOCTET * buffer;      /* start of buffer */
    boolean start_of_file;    /* have we gotten any data yet? */
} my_source_mgr;

typedef my_source_mgr * my_src_ptr;

}

#define INPUT_BUF_SIZE  4096    /* choose an efficiently fread'able size */


/*
 * Initialize source --- called by jpeg_read_header
 * before any data is actually read.
 */

METHODDEF(void)
my_init_source (j_decompress_ptr cinfo)
{
    my_src_ptr src = (my_src_ptr) cinfo->src;

    /* We reset the empty-input-file flag for each image,
     * but we don't clear the input buffer.
     * This is correct behavior for reading a series of images from one source.
     */
    src->start_of_file = TRUE;
}


/*
 * Fill the input buffer --- called whenever buffer is emptied.
 *
 * In typical applications, this should read fresh data into the buffer
 * (ignoring the current state of next_input_byte & bytes_in_buffer),
 * reset the pointer & count to the start of the buffer, and return TRUE
 * indicating that the buffer has been reloaded.  It is not necessary to
 * fill the buffer entirely, only to obtain at least one more byte.
 *
 * There is no such thing as an EOF return.  If the end of the file has been
 * reached, the routine has a choice of ERREXIT() or inserting fake data into
 * the buffer.  In most cases, generating a warning message and inserting a
 * fake EOI marker is the best course of action --- this will allow the
 * decompressor to output however much of the image is there.  However,
 * the resulting error message is misleading if the real problem is an empty
 * input file, so we handle that case specially.
 *
 * In applications that need to be able to suspend compression due to input
 * not being available yet, a FALSE return indicates that no more data can be
 * obtained right now, but more may be forthcoming later.  In this situation,
 * the decompressor will return to its caller (with an indication of the
 * number of scanlines it has read, if any).  The application should resume
 * decompression after it has loaded more data into the input buffer.  Note
 * that there are substantial restrictions on the use of suspension --- see
 * the documentation.
 *
 * When suspending, the decompressor will back up to a convenient restart point
 * (typically the start of the current MCU). next_input_byte & bytes_in_buffer
 * indicate where the restart point will be if the current call returns FALSE.
 * Data beyond this point must be rescanned after resumption, so move it to
 * the front of the buffer rather than discarding it.
 */

METHODDEF(boolean)
my_fill_input_buffer (j_decompress_ptr cinfo)
{
    my_src_ptr src = (my_src_ptr) cinfo->src;
    size_t nbytes;

    nbytes = JFREAD(src->infile, src->buffer, INPUT_BUF_SIZE);

    if (nbytes == 0) {
        if (src->start_of_file) { /* Treat empty input file as fatal error */
            ERREXIT(cinfo, JERR_INPUT_EMPTY);
        }

        WARNMS(cinfo, JWRN_JPEG_EOF);
        /* Insert a fake EOI marker */
        src->buffer[0] = (JOCTET) 0xFF;
        src->buffer[1] = (JOCTET) JPEG_EOI;
        nbytes = 2;
    }

    if (src->start_of_file) {
        src->buffer[0] = (JOCTET) 0xFF;
    }

    src->pub.next_input_byte = src->buffer;
    src->pub.bytes_in_buffer = nbytes;
    src->start_of_file = FALSE;

    return TRUE;
}


/*
 * Skip data --- used to skip over a potentially large amount of
 * uninteresting data (such as an APPn marker).
 *
 * Writers of suspendable-input applications must note that skip_input_data
 * is not granted the right to give a suspension return.  If the skip extends
 * beyond the data currently in the buffer, the buffer can be marked empty so
 * that the next read will cause a fill_input_buffer call that can suspend.
 * Arranging for additional bytes to be discarded before reloading the input
 * buffer is the application writer's problem.
 */

METHODDEF(void)
my_skip_input_data (j_decompress_ptr cinfo, long num_bytes)
{
    my_src_ptr src = (my_src_ptr) cinfo->src;

    /* Just a dumb implementation for now.  Could use fseek() except
     * it doesn't work on pipes.  Not clear that being smart is worth
     * any trouble anyway --- large skips are infrequent.
     */
    if (num_bytes > 0) {
        while (num_bytes > (long) src->pub.bytes_in_buffer) {
            num_bytes -= (long) src->pub.bytes_in_buffer;
            (void) my_fill_input_buffer(cinfo);
            /* note we assume that fill_input_buffer will never return FALSE,
             * so suspension need not be handled.
             */
        }

        src->pub.next_input_byte += (size_t) num_bytes;
        src->pub.bytes_in_buffer -= (size_t) num_bytes;
    }
}


/*
 * An additional method that can be provided by data source modules is the
 * resync_to_restart method for error recovery in the presence of RST markers.
 * For the moment, this source module just uses the default resync method
 * provided by the JPEG library.  That method assumes that no backtracking
 * is possible.
 */


/*
 * Terminate source --- called by jpeg_finish_decompress
 * after all data has been read.  Often a no-op.
 *
 * NB: *not* called by jpeg_abort or jpeg_destroy; surrounding
 * application must deal with any cleanup that should happen even
 * for error exit.
 */

METHODDEF(void)
my_term_source (j_decompress_ptr cinfo)
{
    /* no work necessary here */
}


/*
 * Prepare for input from a stdio stream.
 * The caller must have already opened the stream, and is responsible
 * for closing it after finishing decompression.
 */

GLOBAL(void)
my_jpeg_stdio_src (j_decompress_ptr cinfo, FILE * infile)
{
    my_src_ptr src;

    /* The source object and input buffer are made permanent so that a series
     * of JPEG images can be read from the same file by calling jpeg_stdio_src
     * only before the first one.  (If we discarded the buffer at the end of
     * one image, we'd likely lose the start of the next one.)
     * This makes it unsafe to use this manager and a different source
     * manager serially with the same JPEG object.  Caveat programmer.
     */
    if (cinfo->src == nullptr) { /* first time for this JPEG object? */
        cinfo->src = (struct jpeg_source_mgr *)
                     (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
                             sizeof(my_source_mgr));
        src = (my_src_ptr) cinfo->src;
        src->buffer = (JOCTET *)
                      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
                              INPUT_BUF_SIZE * sizeof(JOCTET));
    }

    src = (my_src_ptr) cinfo->src;
    src->pub.init_source = my_init_source;
    src->pub.fill_input_buffer = my_fill_input_buffer;
    src->pub.skip_input_data = my_skip_input_data;
    src->pub.resync_to_restart = jpeg_resync_to_restart; /* use default method */
    src->pub.term_source = my_term_source;
    src->infile = infile;
    src->pub.bytes_in_buffer = 0; /* forces fill_input_buffer on first read */
    src->pub.next_input_byte = nullptr; /* until buffer loaded */
}

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
    /* Always display the message */
    (*cinfo->err->output_message) (cinfo);

    /* Let the memory manager delete any temp files before we die */
    //jpeg_destroy(cinfo);

    j_decompress_ptr dinfo = (j_decompress_ptr)cinfo;
//  longjmp (((rt_jpeg_error_mgr*)(dinfo->src))->error_jmp_buf, 1);
    longjmp ((reinterpret_cast<rt_jpeg_error_mgr*>(dinfo->src))  ->error_jmp_buf, 1);
}


//const char * const jpeg_std_message_table[] = {
//#include "jerror.h"
//  NULL
//};
extern const char * const jpeg_std_message_table[];


/*
 * Actual output of an error or trace message.
 * Applications may override this method to send JPEG messages somewhere
 * other than stderr.
 *
 * On Windows, printing to stderr is generally completely useless,
 * so we provide optional code to produce an error-dialog popup.
 * Most Windows applications will still prefer to override this routine,
 * but if they don't, it'll do something at least marginally useful.
 *
 * NOTE: to use the library in an environment that doesn't support the
 * C stdio library, you may have to delete the call to fprintf() entirely,
 * not just not use this routine.
 */

METHODDEF(void)
output_message (j_common_ptr cinfo)
{
    char buffer[JMSG_LENGTH_MAX];

    /* Create the message */
    (*cinfo->err->format_message) (cinfo, buffer);

#ifdef USE_WINDOWS_MESSAGEBOX
    /* Display it in a message dialog box */
    MessageBox(GetActiveWindow(), buffer, "JPEG Library Error",
               MB_OK | MB_ICONERROR);
#else
    /* Send it to stderr, adding a newline */
    fprintf(stderr, "%s\n", buffer);
#endif
}


/*
 * Decide whether to emit a trace or warning message.
 * msg_level is one of:
 *   -1: recoverable corrupt-data warning, may want to abort.
 *    0: important advisory messages (always display to user).
 *    1: first level of tracing detail.
 *    2,3,...: successively more detailed tracing messages.
 * An application might override this method if it wanted to abort on warnings
 * or change the policy about which messages to display.
 */

METHODDEF(void)
emit_message (j_common_ptr cinfo, int msg_level)
{
    struct jpeg_error_mgr * err = cinfo->err;

    if (msg_level < 0) {
        /* It's a warning message.  Since corrupt files may generate many warnings,
         * the policy implemented here is to show only the first warning,
         * unless trace_level >= 3.
         */
        if (err->num_warnings == 0 || err->trace_level >= 3) {
            (*err->output_message) (cinfo);
        }

        /* Always count warnings in num_warnings. */
        err->num_warnings++;
    } else {
        /* It's a trace message.  Show it if trace_level >= msg_level. */
        if (err->trace_level >= msg_level) {
            (*err->output_message) (cinfo);
        }
    }
}


/*
 * Format a message string for the most recent JPEG error or message.
 * The message is stored into buffer, which should be at least JMSG_LENGTH_MAX
 * characters.  Note that no '\n' character is added to the string.
 * Few applications should need to override this method.
 */

METHODDEF(void)
format_message (j_common_ptr cinfo, char * buffer)
{
    struct jpeg_error_mgr * err = cinfo->err;
    int msg_code = err->msg_code;
    const char * msgtext = nullptr;
    const char * msgptr;
    char ch;
    boolean isstring;

    /* Look up message string in proper table */
    if (msg_code > 0 && msg_code <= err->last_jpeg_message) {
        msgtext = err->jpeg_message_table[msg_code];
    } else if (err->addon_message_table != nullptr &&
               msg_code >= err->first_addon_message &&
               msg_code <= err->last_addon_message) {
        msgtext = err->addon_message_table[msg_code - err->first_addon_message];
    }

    /* Defend against bogus message number */
    if (msgtext == nullptr) {
        err->msg_parm.i[0] = msg_code;
        msgtext = err->jpeg_message_table[0];
    }

    /* Check for string parameter, as indicated by %s in the message text */
    isstring = FALSE;
    msgptr = msgtext;

    while ((ch = *msgptr++) != '\0') {
        if (ch == '%') {
            if (*msgptr == 's') {
                isstring = TRUE;
            }

            break;
        }
    }

    /* Format the message into the passed buffer */
    if (isstring) {
        sprintf(buffer, msgtext, err->msg_parm.s);
    } else
        sprintf(buffer, msgtext,
                err->msg_parm.i[0], err->msg_parm.i[1],
                err->msg_parm.i[2], err->msg_parm.i[3],
                err->msg_parm.i[4], err->msg_parm.i[5],
                err->msg_parm.i[6], err->msg_parm.i[7]);
}


/*
 * Reset error state variables at start of a new image.
 * This is called during compression startup to reset trace/error
 * processing to default state, without losing any application-specific
 * method pointers.  An application might possibly want to override
 * this method if it has additional error processing state.
 */

METHODDEF(void)
reset_error_mgr (j_common_ptr cinfo)
{
    cinfo->err->num_warnings = 0;
    /* trace_level is not reset since it is an application-supplied parameter */
    cinfo->err->msg_code = 0; /* may be useful as a flag for "no error" */
}


GLOBAL(struct jpeg_error_mgr *)
my_jpeg_std_error (struct jpeg_error_mgr * err)
{

    err->error_exit = my_error_exit;
    err->emit_message = emit_message;
    err->output_message = output_message;
    err->format_message = format_message;
    err->reset_error_mgr = reset_error_mgr;

    err->trace_level = 0;     /* default = no tracing */
    err->num_warnings = 0;    /* no warnings emitted yet */
    err->msg_code = 0;        /* may be useful as a flag for "no error" */

    /* Initialize message table pointers */
    err->jpeg_message_table = jpeg_std_message_table;
    err->last_jpeg_message = (int) JMSG_LASTMSGCODE - 1;

    err->addon_message_table = nullptr;
    err->first_addon_message = 0; /* for safety */
    err->last_addon_message = 0;

    return err;
}
