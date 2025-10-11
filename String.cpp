/*
 Author: Fabio Cannizzo
 NUS FE5226 Example Program

 A dumb String class: just an academic exercise!
*/

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;

class String
{
    char* m_data = nullptr;

    void reset(const char* s)
    {
        if (s != m_data) {
            // release previous memory
            delete[] m_data;
            // allocate new memory: add 1 because of the null terminating character at the end
            m_data = new char[strlen(s) + 1];
            // copy from s
            strcpy(m_data, s);
        }
    }

    friend  ostream& operator<<(ostream& os, const String& s);

public:
    //constructor
    String() { reset(""); }  // points to a string of null size (i.e. just one null character {0}
    String(const char* s) { reset(s); }
    String(const String& s) { reset(s.m_data); }
    //String(const String& s) : String(s.m_data) {}

    // destructor
    ~String() { delete[] m_data; }

    // lenght of the string
    size_t length() const { return strlen(m_data); }

    const char* c_str() const
    {
        return m_data;
    }

    // first overload += operator
    String& operator+=(const char* rhs)
    {
        auto n2 = strlen(rhs);
        if (n2 > 0) {
            auto n1 = length();
            char* tmp = new char[n1 + n2 + 1];  // allocate new memory storage
            strncpy(tmp, m_data, n1);
            strcpy(tmp + n1, rhs);
            delete[] m_data;               // delete old memory storage
            m_data = tmp;                   // point to the new memory storage
        }
        return *this;
    }

    // second overload += operator
    String& operator+=(const String& rhs)
    {
        return *this += rhs.m_data;  // dispatch to the other overload operator +=, which takes an argument const char*
    }

    // first overload operator=
    // makes a copy of the string pointed by rhs
    String& operator=(const char* rhs)
    {
        reset(rhs);
        return *this;
    }

    // second overload operator=
    // makes a deep copy of the string contained in rhs
    String& operator=(const String& rhs)
    {
        return *this = rhs.m_data;   // dispatch to the other overload operator =, which takes an argument const char*
    }
};

// streaming operator
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

