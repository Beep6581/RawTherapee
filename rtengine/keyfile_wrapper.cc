#include "keyfile_wrapper.h"
#include <string>
#include <sstream>

#include <iostream>

namespace rtengine {
    
    const Glib::ustring KeyFileWrapper::INT64_PREFIX("int64");
    
    const Glib::ustring KeyFileWrapper::UINT64_PREFIX("uint64");
    
    const Glib::ustring KeyFileWrapper::BOOL_LIST_PREFIX("bools");
    
    const Glib::ustring KeyFileWrapper::INT_LIST_PREFIX("ints");
    
    const Glib::ustring KeyFileWrapper::DOUBLE_LIST_PREFIX("doubles");
    
    const Glib::ustring KeyFileWrapper::STRING_LIST_PREFIX("strings");
    
    KeyFileWrapper::KeyFileWrapper() {
        // nothing to do here yet
    }
    
    Glib::KeyFile& KeyFileWrapper::getKeyFile() {
        return keyFile;
    }
    
    /**
     * Extracts the string between the first opening and the last closing parenthesis.
     * Example: `extract_from_parentheses("data(Hello world!)")` should return
     * "Hello world!".
     */
    bool KeyFileWrapper::extract_from_parentheses(
            const std::string& s,
            std::string& target,
            std::string open_bracket,
            std::string close_bracket
    ) {
        size_t pos_open = s.find_first_of(open_bracket);
        size_t pos_close = s.find_last_of(close_bracket);
        if (pos_open == std::string::npos || pos_close == std::string::npos) {
            return false;
        }
        target = s.substr(pos_open + 1, pos_close - pos_open - 1);
        return true;
    }
    
    Glib::ustring KeyFileWrapper::read_string_list(const Glib::ustring& in) {
        std::istringstream string_list_stream(in);
        std::ostringstream out;
        
        std::string extracted_item;
        for (std::string item; std::getline(string_list_stream, item, ';');) {
            // TODO: Respect escape sequences like "\;"! (Not needed in current version of RT)
            extract_from_parentheses(item, extracted_item, "\"", "\"");
            out << extracted_item << ';';
        }
        
        return out.str();
    }
    
    bool KeyFileWrapper::load_from_data(const Glib::ustring& data) {
        std::istringstream data_in(data);
        std::ostringstream data_out;
        
        /* Prepare data: Omit everything that's not standard key file syntax. */
        std::string key, value_in, value_out;
        for (std::string line; std::getline(data_in, line);) {
            if (line.find_first_not_of("\t\n\v\f\r") == std::string::npos) {
                continue;
            }
            
            std::istringstream line_stream(line);
            if (!std::getline(line_stream, key, '=') ||
                !std::getline(line_stream, value_in, '=') ||
                (line.size() > 0 && line[0] == '#')
            ) {
                // Line cannot be split into key/value or line is a comment.
                data_out << line << std::endl;
                continue;
            }
            
            // Whitespace around '=' should be ignored.
            // Remove leading/trailing whitespace for easier comparison.
            value_in = g_strstrip(&value_in[0]);
            
            if (value_in.size() > 0 && value_in[0] == '"' &&
                extract_from_parentheses(value_in, value_out, "\"", "\"")
            ) {
                // value_in is formatted as string
                data_out << key << '=' << value_out << std::endl;
                continue;
            }
            
            if (extract_from_parentheses(value_in, value_out)) {
                // value_in is of the form "data_type(content)"
                if (!value_in.compare(0, STRING_LIST_PREFIX.size(), STRING_LIST_PREFIX)) {
                    // content in parentheses is a list of strings, strings still need to be extracted
                    data_out << key << '=' << read_string_list(value_out) << std::endl;
                    continue;
                }
                data_out << key << '=' << value_out << std::endl;
                continue;
            }
            
            // value_in is a regular int, boolean, or a legacy value (e.g. string without quotation marks).
            data_out << key << '=' << value_in << std::endl;
        }
        
        /* Now data_out should look exactly like any other key file content. */
        return keyFile.load_from_data(data_out.str());
    }
    
    Glib::ustring KeyFileWrapper::to_data() {
        return keyFile.to_data();
    }
    
    void KeyFileWrapper::set_boolean(const Glib::ustring& group_name, const Glib::ustring& key, bool value) {
        /* Booleans are stored unchanged. */
        keyFile.set_boolean(group_name, key, value);
    }
    
    void KeyFileWrapper::set_double(const Glib::ustring& group_name, const Glib::ustring& key, double value) {
        /* Doubles must contain a decimal point. */
        std::ostringstream ss;
        ss << std::showpoint << value;
        keyFile.set_string(group_name, key, ss.str());
    }
    
    void KeyFileWrapper::set_integer(const Glib::ustring& group_name, const Glib::ustring& key, int value) {
        /* Integers are stored unchanged. */
        keyFile.set_integer(group_name, key, value);
    }
    
    void KeyFileWrapper::set_int64(const Glib::ustring& group_name, const Glib::ustring& key, gint64 value) {
        /* 64-bit integers must be enclosed with "int64()". */
        std::ostringstream ss;
        ss << INT64_PREFIX << '(' << value << ')';
        keyFile.set_string(group_name, key, ss.str());
    }
    
    void KeyFileWrapper::set_uint64(const Glib::ustring& group_name, const Glib::ustring& key, guint64 value) {
        /* 64-bit integers must be enclosed with "int64()". */
        std::ostringstream ss;
        ss << UINT64_PREFIX << '(' << value << ')';
        keyFile.set_string(group_name, key, ss.str());
    }
    
    void KeyFileWrapper::set_string(const Glib::ustring& group_name, const Glib::ustring& key, const Glib::ustring& string) {
        /* Strings must be surrounded with quotation marks. */
        std::ostringstream ss;
        ss << "\"" << string << "\"";
        keyFile.set_string(group_name, key, ss.str());
    }
    
    void KeyFileWrapper::set_boolean_list(const Glib::ustring& group_name, const Glib::ustring& key, const std::vector<bool>&  list) {
        std::ostringstream ss;
        ss << BOOL_LIST_PREFIX << '(';
        ss << std::boolalpha;
        for (const bool& b : list) {
            ss << b << ";";
        }
        ss << ')';
        keyFile.set_string(group_name, key, ss.str());
    }
    
    void KeyFileWrapper::set_double_list(const Glib::ustring& group_name, const Glib::ustring& key, const std::vector<double>&  list) {
        std::ostringstream ss;
        ss << DOUBLE_LIST_PREFIX << '(';
        ss << std::showpoint;
        for (const double& d : list) {
            ss << d << ";";
        }
        ss << ')';
        keyFile.set_string(group_name, key, ss.str());
    }
    
    void KeyFileWrapper::set_integer_list(const Glib::ustring& group_name, const Glib::ustring& key, const std::vector<int>&  list) {
        std::ostringstream ss;
        ss << INT_LIST_PREFIX << '(';
        for (const int& i : list) {
            ss << i << ";";
        }
        ss << ')';
        keyFile.set_string(group_name, key, ss.str());
    }
    
    void KeyFileWrapper::set_string_list(const Glib::ustring& group_name, const Glib::ustring& key, const std::vector<Glib::ustring>&  list) {
        std::ostringstream ss;
        ss << STRING_LIST_PREFIX << '(';
        for (const Glib::ustring& s : list) {
            ss << "\"" << s << "\";";
        }
        ss << ')';
        keyFile.set_string(group_name, key, ss.str());
    }
    
}
