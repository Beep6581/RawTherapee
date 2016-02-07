#ifndef JAGGEDARRAY_H
#define JAGGEDARRAY_H

namespace rtengine
{

// These emulate a jagged array, but use only 2 allocations instead of 1 + W.

template<class T>
inline T** const allocJaggedArray (const int W, const int H, const bool initZero = false)
{
    T** const a = new T*[H];
    a[0] = new T[H * W];

    for (int i = 1; i < H; ++i) {
        a[i] = a[i - 1] + W;
    }

    if (initZero) {
        std::memset(a[0], 0, sizeof(T) * W * H);
    }

    return a;
}

template<class T>
inline void freeJaggedArray (T** const a)
{
    delete [] a[0];
    delete [] a;
}

template<class T>
class JaggedArray
{
public:
    JaggedArray (const int W, const int H, const bool initZero = false)
    {
        a = allocJaggedArray<T> (W, H, initZero);
    }
    ~JaggedArray ()
    {
        if (a) {
            freeJaggedArray<T> (a);
            a = nullptr;
        }
    }

    JaggedArray (const JaggedArray&) = delete;
    JaggedArray& operator= (const JaggedArray&) = delete;

public:
    operator T** const () const
    {
        return a;
    }

private:
    T** a;

};

// Declared but not defined to prevent
// explicitly freeing a JaggedArray<T> implicitly cast to T**.
template<class T>
void freeJaggedArray (JaggedArray<T>&);

} // rtengine

#endif // JAGGEDARRAY_H
