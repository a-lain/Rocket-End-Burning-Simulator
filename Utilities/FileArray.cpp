#ifndef FILE_ARRAY_CPP
#define FILE_ARRAY_CPP

#include "FileArray.hpp"
#include "FileOperations.hpp"

#include <iostream>

using namespace Utilities;

template <typename T>
unsigned long FileArray<T>::size() const
{
    return N;
}

template <typename T>
T FileArray<T>::operator[](const unsigned long i) const
{
    fseek(file, indexes[i], SEEK_SET);
    T res;
    read_from_file(res, file);
    return res;
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::operator[](const unsigned long i)
{
    return FileArray<T>::Iterator(this, i);
}

template <typename T>
T FileArray<T>::at(const unsigned long i) const
{
    return (*this)[i];
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::at(const unsigned long i)
{
    return (*this)[i];
}

template <typename T>
T FileArray<T>::front() const
{
    return (*this)[0];
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::front()
{
    return (*this)[0];
}

template <typename T>
T FileArray<T>::back() const
{
    return (*this)[N - 1];
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::back()
{
    return (*this)[N - 1];
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::begin()
{
    return FileArray<T>::Iterator(this, 0);
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::end()
{
    return FileArray<T>::Iterator(this, N);
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::rbegin()
{
    return FileArray<T>::Iterator(this, N - 1);
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::rend()
{
    return FileArray<T>::Iterator(this, -1);
}

template <typename T>
void FileArray<T>::push_back(const T& val)
{
    indexes.push_back(0);
    (*this)[N] = val;
    N++;
}

template <typename T>
void FileArray<T>::pop_back()
{
    N--;
    indexes.pop_back();
}

template <typename T>
void FileArray<T>::clear()
{
    N = 0;
    indexes = std::vector<unsigned long>({0});
}

template <typename T>
FileArray<T>::Iterator::Iterator(FileArray<T>* file_array, const unsigned long i):
    file_array(file_array), i(i)
{

}

template <typename T>
bool FileArray<T>::Iterator::operator==(const Iterator& J)
{
    if (this->file_array != J.file_array)
    {
        return false;
    }
    else if (this->i != J.i)
    {
        return false;
    }
    else
    {
        return true;
    }
}

template <typename T>
bool FileArray<T>::Iterator::operator!=(const Iterator& J)
{
    return !((*this) == J);
}

template <typename T>
FileArray<T>::Iterator::operator T() const
{
    fseek(file_array->file, file_array->indexes[i], SEEK_SET);
    T res;
    read_from_file(res, file_array->file);
    return res;
}

template <typename T>
typename FileArray<T>::Iterator& FileArray<T>::Iterator::operator=(const T& val)
{
    fseek(file_array->file,file_array->indexes[i], SEEK_SET);
    write_to_file(val, file_array->file);
    file_array->indexes[i+1] = ftell(file_array->file);
    return *this;
}

template <typename T>
typename FileArray<T>::Iterator& FileArray<T>::Iterator::operator++()
{
    i++;
    return *this;
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::Iterator::operator++(int)
{
    return Iterator(file_array, i++);
}


template <typename T>
typename FileArray<T>::Iterator& FileArray<T>::Iterator::operator--()
{
    i--;
    return *this;
}

template <typename T>
typename FileArray<T>::Iterator FileArray<T>::Iterator::operator--(int)
{
    return Iterator(file_array, i--);
}

template <typename T>
typename FileArray<T>::Iterator& FileArray<T>::Iterator::operator+(const unsigned long j)
{
    i += j;
    return *this;
}

template <typename T>
typename FileArray<T>::Iterator& FileArray<T>::Iterator::operator-(const unsigned long j)
{
    i -= j;
    return *this;
}

template <typename T>
FileArray<T>::FileArray()
{
    this->temp = true;
    this->N = 0;
    char buffer[L_tmpnam];
    tmpnam(buffer);
    this->file = fopen(buffer, "wb+");
    this->filename = std::string(buffer);
    this->indexes = std::vector<unsigned long>({0});
}

template <typename T>
FileArray<T>::FileArray(char* file, const bool temp)
{
    this->temp = temp;
    this->filename = std::string(file);

    // We check if the file exists.
    if (FILE* f = fopen(file, "rb"))
    {
        // We read file length and byte position of all objects.
        read_from_file(this->N, f);
        this->indexes = std::vector<unsigned long>(N + 1);
        for (unsigned long i = 0; i < N + 1; i++)
        {
            read_from_file(indexes[i], f);
        }

        // We create a temporary file.
        char tmp_name[L_tmpnam];
        tmpnam(tmp_name);
        FILE* copy = fopen(tmp_name, "wb+");
        T tmp_val;
        for (unsigned long i = 0; i < N; i++)
        {
            read_from_file(tmp_val, f);
            write_to_file(tmp_val, copy);
        }
        fclose(f);
        rewind(copy);

        // We copy everything back to the main file.
        this->file = fopen(file, "wb+");
        for (unsigned long i = 0; i < N; i++)
        {
            read_from_file(tmp_val, copy);
            write_to_file(tmp_val, this->file);
        }
        fclose(copy);
        system((std::string("rm ") + tmp_name).c_str());
    }
    else
    {
        this->file = fopen(file, "wb+");
        this->N = 0u;
        this->indexes = std::vector<unsigned long>({0u});
    }
}

template <typename T>
FileArray<T>::FileArray(const std::string& file, const bool temp)
{
    this->temp = temp;
    this->filename = std::string(file);

    // We check if the file exists.
    if (FILE* f = fopen(file.c_str(), "rb"))
    {
        // We read file length and byte position of all objects.
        read_from_file(this->N, f);
        this->indexes = std::vector<unsigned long>(N + 1);
        for (unsigned long i = 0; i < N + 1; i++)
        {
            read_from_file(indexes[i], f);
        }

        // We create a temporary file.
        char tmp_name[L_tmpnam];
        tmpnam(tmp_name);
        FILE* copy = fopen(tmp_name, "wb+");
        T tmp_val;
        for (unsigned long i = 0; i < N; i++)
        {
            read_from_file(tmp_val, f);
            write_to_file(tmp_val, copy);
        }
        fclose(f);
        rewind(copy);

        // We copy everything back to the main file.
        this->file = fopen(file.c_str(), "wb+");
        for (unsigned long i = 0; i < N; i++)
        {
            read_from_file(tmp_val, copy);
            write_to_file(tmp_val, this->file);
        }
        fclose(copy);
        system((std::string("rm ") + tmp_name).c_str());
    }
    else
    {
        this->file = fopen(file.c_str(), "wb+");
        this->N = 0;
        this->indexes = std::vector<unsigned long>({0});
    }
}

template <typename T>
FileArray<T>::FileArray(const std::vector<T>& vector)
{
    this->temp = true;
    this->N = 0u;
    this->indexes = std::vector<unsigned long>({0});
    char buffer[L_tmpnam];
    tmpnam(buffer);
    this->file = fopen(buffer, "wb+");
    this->filename = std::string(buffer);

    for (unsigned long j = 0; j < N; j++)
    {
        push_back(vector[j]);
    }
}

template <typename T>
FileArray<T>::FileArray(char* file, const std::vector<T>& vector, const bool temp)
{
    this->temp = temp;
    this->N = 0;
    this->indexes = std::vector<unsigned long>({0});
    this->filename = std::string(file);        
    this->file = fopen(file, "wb+");

    for (unsigned long j = 0; j < N; j++)
    {
        push_back(vector[j]);
    }
}

template <typename T>
FileArray<T>::FileArray(const std::string& file, const std::vector<T>& vector, const bool temp)
{
    this->temp = temp;
    this->N = 0;
    this->indexes = std::vector<unsigned long>({0});
    this->filename = std::string(file);        
    this->file = fopen(file.c_str(), "wb+");

    for (unsigned long j = 0; j < N; j++)
    {
        push_back(vector[j]);
    }
}

template <typename T>
FileArray<T>::FileArray(T* array, unsigned int N)
{
    this->temp = true;
    this->N = 0;
    this->indexes = std::vector<unsigned long>({0});
    char buffer[L_tmpnam];
    tmpnam(buffer);
    this->file = fopen(buffer, "wb+");
    this->filename = std::string(buffer);

    for (unsigned long j = 0; j < N; j++)
    {
        push_back(array[j]);
    }
}

template <typename T>
FileArray<T>::FileArray(char* file, T* array, unsigned int N, const bool temp)
{
    this->temp = temp;
    this->N = 0;
    this->indexes = std::vector<unsigned long>({0});
    this->filename = std::string(file);        
    this->file = fopen(file, "wb+");

    for (unsigned long j = 0; j < N; j++)
    {
        push_back(array[j]);
    }
}

template <typename T>
FileArray<T>::FileArray(const std::string& file, T* array, unsigned int N, const bool temp)
{
    this->temp = temp;
    this->N = 0;
    this->indexes = std::vector<unsigned long>({0});
    this->filename = std::string(file);        
    this->file = fopen(file.c_str(), "wb+");

    for (unsigned long j = 0; j < N; j++)
    {
        push_back(array[j]);
    }
}

template <typename T>
FileArray<T>::FileArray(const std::initializer_list<T>& list)
{
    this->temp = true;
    this->N = 0;
    this->indexes = std::vector<unsigned long>({0});
    char buffer[L_tmpnam];
    tmpnam(buffer);
    this->file = fopen(buffer, "wb+");
    this->filename = std::string(buffer);

    auto it = list.begin();
    for (unsigned long j = 0; j < N; j++)
    {
        push_back(*it);
        ++it;
    }
}

template <typename T>
FileArray<T>::FileArray(char* file, const std::initializer_list<T>& list, const bool temp)
{
    this->temp = temp;
    this->N = 0;
    this->indexes = std::vector<unsigned long>({0});
    this->filename = std::string(file);
    this->file = fopen(file, "wb+");

    auto it = list.begin();
    for (unsigned long j = 0; j < N; j++)
    {
        push_back(*it);
        ++it;
    }
}

template <typename T>
FileArray<T>::FileArray(const std::string& file, const std::initializer_list<T>& list, const bool temp)
{
    this->temp = temp;
    this->N = 0;
    this->indexes = std::vector<unsigned long>({0});
    this->filename = std::string(file);
    this->file = fopen(file.c_str(), "wb+");

    auto it = list.begin();
    for (unsigned long j = 0; j < N; j++)
    {
        push_back(*it);
        ++it;
    }
}

template <typename T>
FileArray<T>::~FileArray()
{
    if (temp || N == 0)
    {
        fclose(file);
        system((std::string("rm ") + filename).c_str());
    }
    else
    {
        // We create a temporal file.
        char tmp_name[L_tmpnam];
        tmpnam(tmp_name);
        FILE* copy = fopen(tmp_name, "wb+");
        for (unsigned long i = 0; i < N; i++)
        {
            write_to_file((T)(*this)[i], copy);
        }
        rewind(copy);

        // We close and reopen the current file to change to saving state.
        fclose(file);
        file = fopen(filename.c_str(), "wb");
        write_to_file(N, file);
        for (unsigned long i = 0; i < N + 1; i++)
        {
            write_to_file(indexes[i], file);
        }
        T tmp_val;
        for (unsigned long i = 0; i < N; i++)
        {
            read_from_file(tmp_val, copy);
            write_to_file(tmp_val, file);
        }
        fclose(copy);
        system((std::string("rm ") + tmp_name).c_str());
        fclose(file);
    }
}

// TEMPLATE INSTANTIATION
namespace Utilities
{
    template class FileArray<short int>;
    template class FileArray<unsigned short int>;
    template class FileArray<int>;
    template class FileArray<unsigned int>;
    template class FileArray<long>;
    template class FileArray<unsigned long>;
    template class FileArray<long long>;
    template class FileArray<unsigned long long>;
    template class FileArray<signed char>;
    template class FileArray<unsigned char>;
    template class FileArray<char>;
    template class FileArray<wchar_t>;
    template class FileArray<float>;
    template class FileArray<double>;
    template class FileArray<long double>;
}

#endif // FILE_ARRAY_CPP
