#ifndef TESTCASE_ITERATOR_H
#define TESTCASE_ITERATOR_H

#include <iostream>
#include <istream>

template <typename TestCase>
class TestCaseIterator
{
public:
    TestCaseIterator() = delete;
    ~TestCaseIterator() {}
    TestCaseIterator(std::istream* const input);
    TestCaseIterator(const TestCaseIterator&) = delete;
    TestCaseIterator(TestCaseIterator&& original);
    bool operator!=(const TestCaseIterator&);
    TestCase& operator*();
    TestCase& operator++();
    TestCase& operator++(int);

    std::istream* input_stream;
    TestCase current;
    int case_count;
    inline void accumulate_cases() { case_count++; }
    virtual TestCase fetch_next() = 0;
};

template <typename TestCase>
TestCaseIterator<TestCase>::TestCaseIterator(std::istream* const input)
    : input_stream(input)
    , current()
    , case_count(0)
{}

template <typename TestCase>
TestCaseIterator<TestCase>::TestCaseIterator(TestCaseIterator&& original)
    : input_stream(original.input_stream)
    , current(std::move(original.current))
{}

template <typename TestCase>
bool TestCaseIterator<TestCase>::operator!=(const TestCaseIterator& rhs)
{
    return input_stream != rhs.input_stream;
}

template <typename TestCase>
TestCase& TestCaseIterator<TestCase>::operator*()
{
    return current;
}

template <typename TestCase>
TestCase& TestCaseIterator<TestCase>::operator++()
{
    current = fetch_next();
    return current;
}

template <typename TestCase>
TestCase& TestCaseIterator<TestCase>::operator++(int)
{
    current = fetch_next();
    return current;
}

#endif  // TESTCASE_ITERATOR_H