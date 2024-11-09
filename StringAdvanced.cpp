/*
Author: Fabio Cannizzo
NUS FE5226 Example Program

A dumb String class

The difference with String.cpp is that for convenience and speed of access
we store information about the length of the string and the capacity of the container
New memory is allocated only if the capacity of the container is not large enough.

*/

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <cstring>

using namespace std;

class String
{
    char* m_data = nullptr;
    size_t m_len = 0;       // length of contained string, exclusive of null terminating character
    size_t m_capacity = 0;  // actual size of the container, inclusive of null terminating character

    void increaseCapacityIfNeeded(size_t newCapacity, bool copyExistent)
    {
        if (m_capacity < newCapacity) {
            // allocate new memory
            char* tmp = new char[newCapacity];
            // if requested, copy existing string to the new block
            if (copyExistent)
                std::copy(m_data, m_data + m_len, tmp);  // do no need to copy the null terminating character
            // release previous memory
            delete[] m_data;
            m_data = tmp;
            m_capacity = newCapacity;
        }
    }

    void replace(const char* s, size_t nNew)
    {
        increaseCapacityIfNeeded(nNew + 1, false);
        strcpy(m_data, s);
        m_len = nNew;
    }

    void append(const char* rhs, size_t n2)
    {
        if (n2) {
            increaseCapacityIfNeeded(m_len + n2 + 1, true);
            strcpy(m_data + m_len, rhs);
            m_len += n2;
        }
    }

    friend  ostream& operator<<(ostream& os, const String& s);

public:
    String() { replace("", 0); }
    String(const char* s) { replace(s, strlen(s)); }
    String(const String& s) { replace(s.m_data, s.length()); }

    ~String() { delete[] m_data; }

    // lenght of the string
    size_t length() const { return m_len; }
    size_t capacity() const { return m_capacity; }

    const char* c_str() const
    {
        return m_data;
    }

    // first overload += operator
    String& operator+=(const char* rhs)
    {
        append(rhs, strlen(rhs));
        return *this;
    }

    // second overload += operator
    String& operator+=(const String& rhs)
    {
        append(rhs.m_data, rhs.length());
        return *this;
    }

    // first overload operator=
    String& operator=(const char* rhs)
    {
        replace(rhs, strlen(rhs));
        return *this;
    }

    // second overload operator=
    String& operator=(const String& rhs)
    {
        replace(rhs.m_data, rhs.length());
        return *this;
    }
};

ostream& operator<<(ostream& os, const String& s)
{
    os << s.m_data;
    return os;
}

int main()
{
    String s1("Hello");
    String s2("World");

    cout << s1 << endl;

    s1 += " ";
    cout << s1 << endl;

    s1 += s2;
    cout << s1 << endl;

    String s3;
    s3 = s1;
    cout << s3 << endl;

    return 0;
}
