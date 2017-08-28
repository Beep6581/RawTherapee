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
#ifndef _KEYFILEWRAPPER_
#define _KEYFILEWRAPPER_

#include <glib.h>
#include <glibmm.h>
#include <vector>

namespace rtengine {

    /**
     * This is a wrapper class for the GNOME glibmm Glib::KeyFile class
     * that reads and outputs strongly typed key files.
     * 
     * Key files saved with the original Glib::KeyFile class do not allow foreign
     * applications to determine the data type of a given key-value pair.
     * For example, saving an integer with value 42 would yield the same result
     * as saving the string "42".
     * This wrapper class makes the following changes to the syntax of key files:
     * 
     * @code
     * # Strings must be enclosed in quotation marks.
     * # (Inner quotation marks need not be escaped, for now.)
     * StringKey="StringValue"
     * # Regular integers are unchanged.
     * IntKey=42
     * # Booleans are unchanged.
     * BoolKey=true
     * # Double values must contain a decimal point. The decimal point need not be
     * # followed by another digit.
     * DoubleKey=42.
     * # 64-bit integers are marked as follows:
     * Int64Key=int64(42)
     * # Similarly, all remaining data types are clearly marked:
     * Uint64Key=uint64(42)
     * BoolListKey=bools(true;false;true;)
     * IntListKey=ints(4;2;42;)
     * DoubleListKey=doubles(4.0;2.0;42.0;)
     * StringListKey=strings("a";"b";)
     * @endcode
     */
    class KeyFileWrapper {
    public:
        static const Glib::ustring
                INT64_PREFIX,
                UINT64_PREFIX,
                BOOL_LIST_PREFIX,
                INT_LIST_PREFIX,
                DOUBLE_LIST_PREFIX,
                STRING_LIST_PREFIX;
        
        /**
         * Creates a new, empty KeyFileWrapper object.
         */
        KeyFileWrapper();
        
        KeyFileWrapper(const KeyFileWrapper&) = delete;
        KeyFileWrapper& operator=(const KeyFileWrapper&) = delete;
        
        /** 
         * Provides access to the wrapped key file.
         * Note that, due to the changes to the syntax, only the getters
         * for strings, booleans, doubles, and integers will work on this object
         * as before.
         * @return The wrapped KeyFile
         */
        Glib::KeyFile& getKeyFile();
        
        /** Loads a KeyFile from memory
         * @param data The data to use as a KeyFile
         * @param flags Bitwise combination of the flags to use for the KeyFile
         * @return true if the KeyFile was successfully loaded, false otherwise
         * @throw Glib::KeyFileError
         */
        bool load_from_data(const Glib::ustring& data);
        
        /** Outputs the KeyFile as a string
         * @return A string object holding the contents of KeyFile
         * @throw Glib::KeyFileError
         */
        Glib::ustring to_data();
        
        void set_boolean(const Glib::ustring& group_name, const Glib::ustring& key, bool value);
        
        void set_double(const Glib::ustring& group_name, const Glib::ustring& key, double value);
        
        void set_integer(const Glib::ustring& group_name, const Glib::ustring& key, int value);

        void set_int64(const Glib::ustring& group_name, const Glib::ustring& key, gint64 value);
        
        void set_uint64(const Glib::ustring& group_name, const Glib::ustring& key, guint64 value);
        
        void set_string(const Glib::ustring& group_name, const Glib::ustring& key, const Glib::ustring& string);
        
        void set_boolean_list(const Glib::ustring& group_name, const Glib::ustring& key, const std::vector<bool>&  list);
        
        void set_double_list(const Glib::ustring& group_name, const Glib::ustring& key, const std::vector<double>&  list);
        
        void set_integer_list(const Glib::ustring& group_name, const Glib::ustring& key, const std::vector<int>&  list);
        
        void set_string_list(const Glib::ustring& group_name, const Glib::ustring& key, const std::vector<Glib::ustring>&  list);
    private:
        Glib::KeyFile keyFile;
        
        /**
         * Returns the substring between the outermost parentheses.
         * Example: `extract_from_parentheses("ints(3;4;5;)")` would return "3;4;5;"
         */
        static bool extract_from_parentheses(
                const std::string& s,
                std::string& target,
                std::string open_bracket="(",
                std::string close_bracket=")"
        );
        
        /**
         * Converts a string list like
         * `"a";"hello";"world";`
         * to a Glib::KeyFile conformant string list like
         * `a;hello;world`.
         */
        static Glib::ustring read_string_list(const Glib::ustring& in);
    };

}

#endif
