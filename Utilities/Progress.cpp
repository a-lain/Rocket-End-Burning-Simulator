#ifndef PROGRESS_CPP
#define PROGRESS_CPP
#include "Progress.hpp"
#include <cmath>
#include <future>
#include <thread>
#include <chrono>
#include <algorithm>
#include <iostream>

using namespace Utilities;

// BEGIN LIBRARY CONFIGURATION
const char* dateFormat = "%c";
// END LIBRARY CONFIGURATION

// BEGIN IMPLEMENTATION

std::string Utilities::time_to_string(const std::chrono::system_clock::time_point date)
{
    static char buffer[100];
    time_t time = std::chrono::system_clock::to_time_t(date);
    strftime(buffer, 100, dateFormat, localtime(&time));
    return std::string(buffer);
}

std::string Utilities::duration_to_string(unsigned int seconds)
{
    unsigned int days = seconds / 86400;
    seconds -= days*86400;
    unsigned int hours = seconds / 3600;
    seconds -= hours*3600;
    unsigned int minutes = seconds / 60;
    seconds -= minutes*60;
    static char buffer[20];
    sprintf(buffer, "%ud %02uh %02um %02us", days, hours, minutes, seconds);
    return std::string(buffer);
}

std::string Utilities::duration_to_string(const double seconds)
{
    return Utilities::duration_to_string((unsigned int) round(seconds));
}

std::string Utilities::duration_to_string(const std::chrono::system_clock::duration time)
{
    return Utilities::duration_to_string((unsigned int) std::chrono::duration_cast<std::chrono::seconds>(time).count());
}

void ProgressBase::start()
{
    status = ProgressStatus::RUNNING;
    t_start = std::chrono::system_clock::now();
}

void ProgressBase::pause()
{
    status = ProgressStatus::PAUSED;
    t_paused = std::chrono::system_clock::now();
}

void ProgressBase::resume()
{
    status = ProgressStatus::RUNNING;
    inactive_time += std::chrono::system_clock::now() - t_paused;
}

void ProgressBase::finish()
{
    status = ProgressStatus::FINISHED;
    t_finish = std::chrono::system_clock::now();
}

std::string ProgressBase::report(unsigned int level) const
{
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::string pre = std::string(level, '\t');
    std::string text = pre + "• " + name + "\n";
    if (status == ProgressStatus::HOLDING)
    {
        text += pre + "  STATUS: HOLDING\n";
    }
    else if (status == ProgressStatus::RUNNING)
    {   
        text += pre + "  STATUS: RUNNING\n";
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Time active elapsed: " + duration_to_string(now - t_start - inactive_time) + "\n";
    }
    else if (status == ProgressStatus::PAUSED)
    {
        text += pre + "  STATUS: PAUSED\n";
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Time active elapsed: " + duration_to_string(now - t_start - inactive_time) + "\n";
    }
    else if (status == ProgressStatus::FINISHED)
    {
        text += pre + "  STATUS: FINISHED\n";
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Finishing time: " + time_to_string(t_finish) + "\n";
        text += pre + "  Duration active time: " + duration_to_string(t_finish - t_start - inactive_time) + "\n";
    }
    for (auto it = children.begin(); it != children.end(); ++it)
    {
        text += (*it->second).report(level + 1);
    }
    return text;
}

ProgressBase::ProgressBase(const std::string& name, const ProgressEstimation estimation)
{
    this->name = name;
    this->status = ProgressStatus::HOLDING;
    this->estimation = estimation;
    this->inactive_time = std::chrono::system_clock::duration(0);
}

ProgressBase::ProgressBase(const ProgressBase& progress)
{
    this->name = progress.name;
    this->estimation = progress.estimation;
    this->t_start = progress.t_start;
    this->t_paused = progress.t_paused;
    this->t_finish = progress.t_finish;
    this->inactive_time = progress.inactive_time;
    this->children = progress.children;
    this->status = progress.status;
}

ProgressBase& ProgressBase::operator=(const ProgressBase& progress)
{
    this->name = progress.name;
    this->estimation = progress.estimation;
    this->t_start = progress.t_start;
    this->t_paused = progress.t_paused;
    this->t_finish = progress.t_finish;
    this->inactive_time = progress.inactive_time;
    this->children = progress.children;
    this->status = progress.status;
    return *this;
}

ProgressBase& ProgressBase::operator[] (std::string name)
{
    return *children[name];
}

void ProgressBase::add_child(ProgressBase& progress)
{
    M.lock();
    children[progress.name] = &progress;
    M.unlock();
}

void ProgressBase::eliminate_child(ProgressBase& progress)
{
    M.lock();
    children.erase(progress.name);
    M.unlock();
}

void ProgressBase::update_to_terminal(unsigned int period)
{
    auto f = [&](unsigned int period)
    {
        unsigned int N_lines_old = 0;
        std::this_thread::sleep_for(std::chrono::milliseconds(period));
        while(status != ProgressStatus::FINISHED)
        {
            // We eliminate unneeded old lines.
            for (unsigned int i = 0; i < N_lines_old; i++)
            {
                std::cout << "\033[A\r\033[2K";
            }

            std::string text = report(0);
            N_lines_old = std::count(text.begin(), text.end(), '\n');
            std::cout << text;

            std::this_thread::sleep_for(std::chrono::milliseconds(period));
        }

        // We show the progress one last time.
        for (unsigned int i = 0; i < N_lines_old; i++)
        {
            std::cout << "\033[A\r\033[2K";
        }
        std::cout << report(0);
    };

    task = std::async(std::launch::async, f, period);
}

template<class T>
Progress<T>::Progress(const std::string& name, ProgressEstimation estimation, const T initial, const T objective) :
    ProgressBase(name, estimation)
{
    this->initial = initial;
    this->current = initial;
    this->objective = objective;
}

template<class T>
Progress<T>::Progress(Progress<T>& progress) : ProgressBase(progress)
{
    this->initial = progress.initial;
    this->current = T(progress.current);
    this->objective = progress.objective;
}

template<class T>
Progress<T>& Progress<T>::operator=(const Progress<T>& progress)
{
    this->ProgressBase::operator=(progress);
    this->initial = progress.initial;
    this->current = T(progress.current);
    this->objective = progress.objective;
    return *this;
}

template<class T>
Progress<T>& Progress<T>::operator=(const T& val)
{
    this->current = val;
    return *this;
}

template <>
Progress<unsigned int>& Progress<unsigned int>::operator+=(const unsigned int& val)
{
    this->current += val;
    return *this;
}

template <>
Progress<unsigned long long>& Progress<unsigned long long>::operator+=(const unsigned long long& val)
{
    this->current += val;
    return *this;
}

template<>
Progress<unsigned int>& Progress<unsigned int>::operator++()
{
    this->current++;
    return *this;
}

template<>
Progress<unsigned long long>& Progress<unsigned long long>::operator++()
{
    this->current++;
    return *this;
}

template<>
Progress<unsigned int>& Progress<unsigned int>::operator-=(const unsigned int& val)
{
    this->current -= val;
    return *this;
}

template<>
Progress<unsigned long long>& Progress<unsigned long long>::operator-=(const unsigned long long& val)
{
    this->current -= val;
    return *this;
}

template<>
Progress<unsigned int>& Progress<unsigned int>::operator--()
{
    this->current--;
    return *this;
}

template<>
Progress<unsigned long long>& Progress<unsigned long long>::operator--()
{
    this->current--;
    return *this;
}

template<class T>
Progress<T>& Progress<T>::operator*=(const T& val)
{
    this->current = this->current * val;
    return *this;
}

template<class T>
Progress<T>& Progress<T>::operator/=(const T& val)
{
    this->current = this->current / val;
    return *this;
}

template<>
Progress<unsigned int>& Progress<unsigned int>::operator%=(const unsigned int& val)
{
    this->current = this->current % val;
    return *this;
}

template<>
Progress<unsigned long long>& Progress<unsigned long long>::operator%=(const unsigned long long& val)
{
    this->current = this->current % val;
    return *this;
}

template<class T>
Progress<T>::operator std::atomic<T>&()
{
    return this->current;
}

template<class T>
Progress<T>::operator T() const
{
    return T(this->current);
}

template<class T>
std::string Progress<T>::report(const unsigned int level) const
{
    return this->ProgressBase::report(level);
}

template<>
std::string Progress<unsigned int>::report(unsigned int level) const
{
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::string pre = std::string(level, '\t');
    std::string text = pre + "• " + name + "\n";
    unsigned int total = objective - initial;
    if (status == ProgressStatus::HOLDING)
    {
        text += pre + "  STATUS: HOLDING\n";
    }
    else if (status == ProgressStatus::RUNNING || status == ProgressStatus::PAUSED)
    {   
        // Delta expressed in seconds
        double delta_t = std::chrono::duration_cast<std::chrono::seconds>(now - t_start - inactive_time).count();
        unsigned int done = current - initial;
        unsigned int remaining = objective - current;
        text += pre + ( (status == ProgressStatus::RUNNING) ? "  STATUS: RUNNING\n" : "  STATUS: PAUSED\n" );
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Time active elapsed: " + duration_to_string(delta_t) + "\n";
        if (delta_t > 0 && total > 0)
        {            
            text += pre + "  Progress: " + std::to_string(done) + "/" + std::to_string(total) +
                "    " + std::to_string( 100 * (double) done / total) + "% complete\n";
            text += pre + "  Speed: " + std::to_string((double) done / delta_t) + " [steps / s]" +
                "    " + std::to_string( 100 * done / total / delta_t) + " [% / s]\n";
            // Estimated time remaining is (steps remaining) / (steps completed) * (Delta t)
            text += pre + "  Estimated active time remaining: " + duration_to_string(
                (double) remaining / done * delta_t ) + "\n";
        }
    }
    else if (status == ProgressStatus::FINISHED)
    {
        text += pre + "  STATUS: FINISHED\n";
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Finishing time: " + time_to_string(t_finish) + "\n";
        text += pre + "  Number of steps completed: " + std::to_string(total) + "\n";
        text += pre + "  Duration active time: " + duration_to_string(t_finish - t_start - inactive_time) + "\n";
    }
    for (auto it = children.begin(); it != children.end(); ++it)
    {
        text += (*it->second).report(level + 1);
    }
    return text;
}

template<>
std::string Progress<unsigned long long int>::report(unsigned int level) const
{
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::string pre = std::string(level, '\t');
    std::string text = pre + "• " + name + "\n";
    unsigned long long int total = objective - initial;
    if (status == ProgressStatus::HOLDING)
    {
        text += pre + "  STATUS: HOLDING\n";
    }
    else if (status == ProgressStatus::RUNNING || status == ProgressStatus::PAUSED)
    {   
        // Delta expressed in seconds
        double delta_t = std::chrono::duration_cast<std::chrono::seconds>(now - t_start - inactive_time).count();
        unsigned long long int done = current - initial;
        unsigned long long int remaining = objective - current;
        text += pre + ( (status == ProgressStatus::RUNNING) ? "  STATUS: RUNNING\n" : "  STATUS: PAUSED\n" );
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Time active elapsed: " + duration_to_string(delta_t) + "\n";
        if (delta_t > 0 && total > 0)
        {            
            text += pre + "  Progress: " + std::to_string(done) + "/" + std::to_string(total) +
                "    " + std::to_string( 100 * (double) done / total) + "% complete\n";
            text += pre + "  Speed: " + std::to_string((double) done / delta_t) + " [steps / s]" +
                "    " + std::to_string( 100 * done / total / delta_t) + " [% / s]\n";
            // Estimated time remaining is (steps remaining) / (steps completed) * (Delta t)
            text += pre + "  Estimated active time remaining: " + duration_to_string(
                (double) remaining / done * delta_t ) + "\n";
        }
    }
    else if (status == ProgressStatus::FINISHED)
    {
        text += pre + "  STATUS: FINISHED\n";
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Finishing time: " + time_to_string(t_finish) + "\n";
        text += pre + "  Number of steps completed: " + std::to_string(total) + "\n";
        text += pre + "  Duration active time: " + duration_to_string(t_finish - t_start - inactive_time) + "\n";
    }
    for (auto it = children.begin(); it != children.end(); ++it)
    {
        text += (*it->second).report(level + 1);
    }
    return text;
}

template<>
std::string Progress<double>::report(unsigned int level) const
{
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::string pre = std::string(level, '\t');
    std::string text = pre + "• " + name + "\n";
    if (status == ProgressStatus::HOLDING)
    {
        text += pre + "  STATUS: HOLDING\n";
    }
    else if (status == ProgressStatus::RUNNING || status == ProgressStatus::PAUSED)
    {   
        double delta_t = std::chrono::duration_cast<std::chrono::seconds>(now - t_start - inactive_time).count();
        text += pre + ( (status == ProgressStatus::RUNNING) ? "  STATUS: RUNNING\n" : "  STATUS: PAUSED\n" );
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Time elapsed: " + duration_to_string(delta_t) + "\n";
        double done, remaining, total;
        if (estimation == ProgressEstimation::LINEAR)
        {
            done = fabs(current - initial);
            remaining = fabs(objective - current);
            total = fabs(objective - initial);
        }
        else // (estimation == LOGARITHMIC)
        {
            done = fabs(log10(current / initial));
            remaining = fabs(log10(objective / current));
            total = fabs(log10(objective / initial));
        }
        if (delta_t > 0 && total > 0)
        {
            text += pre + "  Progress: " + std::to_string(100 * done / total) + "% complete\n";
            text += pre + "  Speed: " + std::to_string( 100 * done / total / delta_t) + " [% / s]\n";
            text += pre + "  Estimated active time remaining: " + duration_to_string(
            remaining /  done  * delta_t ) + "\n";
        }
    }
    else if (status == ProgressStatus::FINISHED)
    {
        text += pre + "  STATUS: FINISHED\n";
        text += pre + "  Starting time: " + time_to_string(t_start) + "\n";
        text += pre + "  Finishing time: " + time_to_string(t_finish) + "\n";
        text += pre + "  Duration active time: " + duration_to_string(t_finish - t_start - inactive_time) + "\n";
    }
    for (auto it = children.begin(); it != children.end(); ++it)
    {
        text += (*it->second).report(level + 1);
    }
    return text;
}


// END IMPLEMENTATION

namespace Utilities
{
    template class Progress<unsigned int>;
    template class Progress<unsigned long long>;
    template class Progress<double>;
}

#endif // PROGRESS_CPP
