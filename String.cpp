/*
 Author: Fabio Cannizzo
 NUS FE5226 Example Program

 A dumb String class
*/

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>

using namespace std;

class String
{
    char *m_data;

    void reset(const char *s)
    {
        if (s != m_data) {
            // release previous memory
            delete [] m_data;
            // allocate new memory: add 1 for the zero at the end
            m_data = new char[strlen(s) + 1];
            // copy from s
            strcpy(m_data, s);
        }
    }

    friend  ostream& operator<<(ostream& os, const String& s);

public:
    String() : m_data(NULL) { reset(""); }
    String(const char *s) : m_data(NULL) { reset(s); }
    String(const String& s) : m_data(NULL) { reset(s.m_data); }

    ~String() { delete[] m_data; }

    // lenght of the string
    unsigned length() const { return strlen(m_data); }

    const char *operator()() const
    {
        return m_data;
    }

    // first overload += operator
    String& operator+=(const char* rhs)
    {
        unsigned n2 = strlen(rhs);
        if (n2>0) {
            unsigned n1 = length();
            char *tmp = new char[n1+n2+1];
            strncpy(tmp, m_data, n1);
            strcpy(tmp+n1, rhs);
            delete [] m_data;
            m_data = tmp;
        }
        return *this;
    }

    // second overload += operator
    String& operator+=(const String& rhs)
    {
        return *this += rhs.m_data;
    }

    // first overload += operator
    String& operator=(const char* rhs)
    {
        reset(rhs);
        return *this;
    }

    // second overload += operator
    String& operator=(const String& rhs)
    {
        return *this = rhs.m_data;
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
