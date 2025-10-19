#ifndef TESTCASE_LOADER_H
#define TESTCASE_LOADER_H

#include <fstream>
#include <iostream>
#include <string>

template <typename Iterator>
class TestCaseLoader
{
private:
    std::fstream m_filestream;
    std::istream* m_datastream;

public:
    TestCaseLoader(const std::string& filename)
    {
        m_filestream.open(filename);
        m_datastream = &m_filestream;
    }
    ~TestCaseLoader() { m_filestream.close(); }
    Iterator begin() { return Iterator{m_datastream}; }
    Iterator end() { return Iterator{nullptr}; }
};

#endif  // TESTCASE_LOADER_H