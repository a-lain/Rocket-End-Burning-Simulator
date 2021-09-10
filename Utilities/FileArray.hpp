#ifndef FILE_ARRAY_HPP
#define FILE_ARRAY_HPP

#include <cstdio>
#include <initializer_list>
#include <string>
#include <vector>

namespace Utilities
{
    template <typename T>
    class FileArray
    {
        public:
        FILE* file;

        /**
         * @brief Returns number of elements in the array.
         * 
         * @return unsigned long 
         */
        unsigned long size() const;

        class Iterator
        {
            public:

            FileArray* file_array;
            unsigned long i;

            explicit Iterator(FileArray* file_array, const unsigned long i = 0);

            bool operator==(const Iterator& J);
            bool operator!=(const Iterator& J);

            operator T() const;

            Iterator& operator=(const T& val);

            Iterator& operator++();
            Iterator operator++(int);
            Iterator& operator--();
            Iterator operator--(int);
            Iterator& operator+(const unsigned long j);
            Iterator& operator-(const unsigned long j);
        };

        T operator[](const unsigned long i) const;
        Iterator operator[](const unsigned long i);
        T at(const unsigned long i) const;
        Iterator at(const unsigned long i);
        T front() const;
        Iterator front();
        T back() const;
        Iterator back();

        Iterator begin();
        Iterator end();
        Iterator rbegin();
        Iterator rend();

        void push_back(const T& val);
        void pop_back();

        // void insert(const unsigned long i, const T& val);
        // void insert(const unsigned long i, const unsigned long n, const T& val);
        // void insert(const unsigned long i, const std::initializer_list<T>& list);
        // void erase(const unsigned long i);
        // void erase(const unsigned long start, unsigned long stop);
        void clear();

        explicit FileArray();
        explicit FileArray(char* file, const bool temp = false);
        explicit FileArray(const std::string& file, const bool temp = false);

        explicit FileArray(const std::vector<T>& vector);
        FileArray(char* file, const std::vector<T>& vector, const bool temp = false);
        FileArray(const std::string& file, const std::vector<T>& vector, const bool temp = false);

        FileArray(T* array, unsigned int N);
        FileArray(char* file, T* array, unsigned int N, const bool temp = false);
        FileArray(const std::string& file, T* array, unsigned int N, const bool temp = false);

        explicit FileArray(const std::initializer_list<T>& list);
        FileArray(char* file, const std::initializer_list<T>& list, const bool temp = false);
        FileArray(const std::string& file, const std::initializer_list<T>& list, const bool temp = false);

        // There is no copy constructor and assigment operator.
        FileArray(const FileArray<T>& f_array) = delete;
        FileArray& operator=(const FileArray<T>& f_array) = delete;
        ~FileArray();

        protected:

        /**
         * @brief Length of the array.
         * 
         */
        unsigned long N;

        std::string filename;

        bool temp;

        std::vector<unsigned long> indexes;
    };
}

#endif // FILE_ARRAY_HPP