#ifndef FILE_OPERATIONS_HPP
#define FILE_OPERATIONS_HPP

#include <complex>
#include <vector>
#include <cstdio>
#include <string>

template <class T>
void inline read_from_file(T& object, FILE* file)
{
    object.read_from_file(file);
}

void inline read_from_file(char& val, FILE* file)
{
    fread(&val, sizeof(char), 1, file);
}

void inline read_from_file(signed char& val, FILE* file)
{
    fread(&val, sizeof(signed char), 1, file);
}

void inline read_from_file(unsigned char& val, FILE* file)
{
    fread(&val, sizeof(unsigned char), 1, file);
}

void inline read_from_file(short int& val, FILE* file)
{
    fread(&val, sizeof(short int), 1, file);
}

void inline read_from_file(unsigned short int& val, FILE* file)
{
    fread(&val, sizeof(unsigned short int), 1, file);
}

void inline read_from_file(int& val, FILE* file)
{
    fread(&val, sizeof(int), 1, file);
}

void inline read_from_file(unsigned int& val, FILE* file)
{
    fread(&val, sizeof(unsigned int), 1, file);
}

void inline read_from_file(long int& val, FILE* file)
{
    fread(&val, sizeof(long int), 1, file);
}

void inline read_from_file(unsigned long int& val, FILE* file)
{
    fread(&val, sizeof(unsigned long int), 1, file);
}

void inline read_from_file(long long int& val, FILE* file)
{
    fread(&val, sizeof(long long int), 1, file);
}

void inline read_from_file(unsigned long long int& val, FILE* file)
{
    fread(&val, sizeof(unsigned long long int), 1, file);
}

void inline read_from_file(float& val, FILE* file)
{
    fread(&val, sizeof(float), 1, file);
}

void inline read_from_file(double& val, FILE* file)
{
    fread(&val, sizeof(double), 1, file);
}

void inline read_from_file(long double& val, FILE* file)
{
    fread(&val, sizeof(long double), 1, file);
}

void inline read_from_file(wchar_t& val, FILE* file)
{
    fread(&val, sizeof(wchar_t), 1, file);
}

template <class K>
void inline read_from_file(std::complex<K>& z, FILE* file)
{
    K x, y;
    read_from_file(x, file);
    read_from_file(y, file);
    z = std::complex<K>(x, y);
}

void inline read_from_file(std::string& val, FILE* file)
{
    char buffer[1024];
    unsigned int i = 0;
    fread(buffer, sizeof(char), 1, file);
    while(buffer[i] != '\0')
    {
        i++;
        fread(buffer + i, sizeof(char), 1, file);
    }
    val = std::string(buffer);
}

template <class T>
void inline read_from_file(std::vector<T>& object, FILE* file)
{
    unsigned int size;
    read_from_file(size, file);
    object = std::vector<T>(size);
    for (unsigned int i = 0; i < size; i++)
    {
        read_from_file(object[i], file);
    }
}

void inline read_from_file(bool& val, FILE* file)
{
    char temp;
    read_from_file(temp, file);
    val = temp;
}


template <class T>
void inline write_to_file(const T& object, FILE* file)
{
    object.write_to_file(file);
}

void inline write_to_file(const signed char& val, FILE* file)
{
    fwrite(&val, sizeof(signed char), 1, file);
}

void inline write_to_file(const char& val, FILE* file)
{
    fwrite(&val, sizeof(char), 1, file);
}

void inline write_to_file(const unsigned char& val, FILE* file)
{
    fwrite(&val, sizeof(unsigned char), 1, file);
}

void inline write_to_file(const short int& val, FILE* file)
{
    fwrite(&val, sizeof(short int), 1, file);
}

void inline write_to_file(const unsigned short int& val, FILE* file)
{
    fwrite(&val, sizeof(unsigned short int), 1, file);
}

void inline write_to_file(const int& val, FILE* file)
{
    fwrite(&val, sizeof(int), 1, file);
}

void inline write_to_file(const unsigned int& val, FILE* file)
{
    fwrite(&val, sizeof(unsigned int), 1, file);
}

void inline write_to_file(const long int& val, FILE* file)
{
    fwrite(&val, sizeof(long int), 1, file);
}

void inline write_to_file(const unsigned long int& val, FILE* file)
{
    fwrite(&val, sizeof(unsigned long int), 1, file);
}

void inline write_to_file(const long long int& val, FILE* file)
{
    fwrite(&val, sizeof(long long int), 1, file);
}

void inline write_to_file(const unsigned long long int& val, FILE* file)
{
    fwrite(&val, sizeof(unsigned long long int), 1, file);
}

void inline write_to_file(const float& val, FILE* file)
{
    fwrite(&val, sizeof(float), 1, file);
}

void inline write_to_file(const double& val, FILE* file)
{
    fwrite(&val, sizeof(double), 1, file);
}

void inline write_to_file(const long double& val, FILE* file)
{
    fwrite(&val, sizeof(long double), 1, file);
}

void inline write_to_file(const wchar_t& val, FILE* file)
{
    fwrite(&val, sizeof(wchar_t), 1, file);
}

template <class K>
void inline write_to_file(const std::complex<K>& z, FILE* file)
{
    K x = z.real();
    K y = z.imag();
    fwrite(&x, sizeof(K), 1, file);
    fwrite(&y, sizeof(K), 1, file);
}

void inline write_to_file(const std::string& val, FILE* file)
{
    fwrite(val.c_str(), sizeof(char), val.size()+1, file);
}

template <class T>
void inline write_to_file(const std::vector<T>& object, FILE* file)
{
    write_to_file((unsigned int)object.size(), file);
    for (unsigned int i = 0; i < object.size(); i++)
    {
        write_to_file(object[i], file);
    }
}

void inline write_to_file(const bool& val, FILE* file)
{
    write_to_file((char) val, file);
}

#endif // FILE_OPERATIONS_HPP