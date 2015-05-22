#ifndef SAFE_KEY_FILE_H_INCLUDED
#define SAFE_KEY_FILE_H_INCLUDED

#include <glibmm.h>
namespace rtengine {
	
class SafeKeyFile : public  Glib::KeyFile
{
	public :
	
#ifdef GLIBMM_EXCEPTIONS_ENABLED
#define SAFE_KEY_FILE_METHOD_CODE(method,method_err) \
			do { try { res = Glib::KeyFile::method; }catch (const Glib::KeyFileError& e) { }	; \
					 return res; }while(0)
#else
#define SAFE_KEY_FILE_METHOD_CODE(method,method_err) \
			do { std::auto_ptr<Glib::Error> error; \
				res = Glib::KeyFile::method_err; \
				if (error.get()){/* TODO */}; \
				return res;} while(0)
#endif //GLIBMM_EXCEPTIONS_ENABLED
#define SAFE_KEY_FILE_METHOD(method,method_err,ret_type) \
		{	ret_type res = (ret_type)0; SAFE_KEY_FILE_METHOD_CODE(method,method_err);}

#define SAFE_KEY_FILE_METHOD_NOINIT(method,method_err,ret_type) \
		{	ret_type res; SAFE_KEY_FILE_METHOD_CODE(method,method_err);}
	
		Glib::ustring to_data()
			SAFE_KEY_FILE_METHOD_NOINIT(to_data(), to_data(error), Glib::ustring);
		
		bool load_from_data(const Glib::ustring& data, Glib::KeyFileFlags flags = Glib::KEY_FILE_NONE)
			SAFE_KEY_FILE_METHOD(load_from_data(data,flags), load_from_data(data,flags,error), bool);
		
		bool load_from_file(const std::string& filename, Glib::KeyFileFlags flags = Glib::KEY_FILE_NONE)		
			SAFE_KEY_FILE_METHOD(load_from_file(filename,flags), load_from_file(filename,flags,error), bool);
		
		bool has_key(const Glib::ustring& group_name, const Glib::ustring& key) const		
			SAFE_KEY_FILE_METHOD(has_key(group_name,key), has_key(group_name,key,error), bool);
			
		bool get_boolean(const Glib::ustring& group_name, const Glib::ustring& key) const
			SAFE_KEY_FILE_METHOD(get_boolean(group_name,key), get_boolean(group_name,key,error), bool);
		
		int get_integer(const Glib::ustring& group_name, const Glib::ustring& key) const
			SAFE_KEY_FILE_METHOD(get_integer(group_name,key), get_integer(group_name,key,error), int);


/*
		double get_double(const Glib::ustring& group_name, const Glib::ustring& key) const
			SAFE_KEY_FILE_METHOD(get_double(group_name,key), get_double(group_name,key,error), double);
			
		typedef std::vector<double> DoubleArrayType;		
	
		DoubleArrayType get_double_list(const Glib::ustring& group_name, const Glib::ustring& key) const
			SAFE_KEY_FILE_METHOD_NOINIT(get_double_list(group_name,key), get_double_list(group_name,key,error), DoubleArrayType);
*/
		typedef std::vector<int> IntArrayType;		
	
		IntArrayType get_integer_list(const Glib::ustring& group_name, const Glib::ustring& key) const
			SAFE_KEY_FILE_METHOD_NOINIT(get_integer_list(group_name,key), get_integer_list(group_name,key,error), IntArrayType);

		Glib::ustring get_string(const Glib::ustring& group_name, const Glib::ustring& key) const
			SAFE_KEY_FILE_METHOD_NOINIT(get_string(group_name,key), get_string(group_name,key,error), Glib::ustring);

		double get_double(const Glib::ustring& group_name, const Glib::ustring& key) const {
		    Glib::ustring temp = get_string( group_name, key);
		    if(!temp.empty()) {
		        double tmpdbl;
                if(sscanf(temp.c_str(), "%lf", &tmpdbl))
                    return tmpdbl;
                else
                    return 0.0;
		    }
            return 0.0;
		}

		typedef std::vector<Glib::ustring> StringArrayType;
			
		StringArrayType get_string_list(const Glib::ustring& group_name, const Glib::ustring& key) const
			SAFE_KEY_FILE_METHOD_NOINIT(get_string_list(group_name,key), get_string_list(group_name,key,error), StringArrayType);

		typedef std::vector<double> DoubleArrayType;		
	
		DoubleArrayType get_double_list(const Glib::ustring& group_name, const Glib::ustring& key) const {
            StringArrayType temp = get_string_list(group_name, key);
            DoubleArrayType tempdouble;
            unsigned int n = temp.size();
            if(n) {
                tempdouble.reserve(n);
          	    for (unsigned int i=0; i<n; i++) {
                    if(!temp[i].empty()) {
                        double tmpdbl;
                        if(sscanf(temp[i].c_str(), "%lf", &tmpdbl))
                            tempdouble.push_back(tmpdbl);
                        else
                            tempdouble.push_back(0.0);
                    } else {
                        tempdouble.push_back(0.0);
                    }
          	    }
            }
            return tempdouble;
		}


		StringArrayType get_keys(const Glib::ustring& group_name) const
			SAFE_KEY_FILE_METHOD_NOINIT(get_keys(group_name), get_keys(group_name,error), StringArrayType);

#undef SAFE_KEY_FILE_METHOD_CODE
#undef SAFE_KEY_FILE_METHOD
#undef SAFE_KEY_FILE_METHOD_NOINIT

};

}

#endif
