#ifndef PROGRESS_HPP
#define PROGRESS_HPP

#include <string>
#include <vector>
#include <atomic>
#include <chrono>
#include <map>
#include <mutex>
#include <future>

namespace Utilities
{
    std::string time_to_string(const std::chrono::system_clock::time_point date);
    std::string duration_to_string(unsigned int seconds);
    std::string duration_to_string(const double seconds);
    std::string duration_to_string(const std::chrono::system_clock::duration time);

    enum class ProgressStatus
    {
        HOLDING,
        RUNNING,
        PAUSED,
        FINISHED,
    };

    enum class ProgressEstimation
    {
        LINEAR,
        LOGARITHMIC,
    };

    class ProgressBase
    {
        protected:
        ProgressEstimation estimation;
        std::string name;
        std::chrono::system_clock::time_point t_start;
        std::chrono::system_clock::time_point t_paused;
        std::chrono::system_clock::time_point t_finish;
        std::chrono::system_clock::duration inactive_time;
        std::map<std::string, ProgressBase*> children;
        std::mutex M;
        std::future<void> task;

        public:
        ProgressStatus status;
        void start();
        void pause();
        void resume();
        void finish();
        virtual std::string report(const unsigned int level) const;

        ProgressBase(const std::string& name = "", const ProgressEstimation estimation = ProgressEstimation::LINEAR);
        ProgressBase(const ProgressBase& progress);
        ProgressBase& operator=(const ProgressBase& progress);

        ProgressBase& operator[] (std::string name);

        void add_child(ProgressBase& progress);
        void eliminate_child(ProgressBase& progress);

        void update_to_terminal(unsigned int period = 250);
    };


    template<class T>
    class Progress: public ProgressBase
    {
        protected:
        T initial;
        T objective;
        std::atomic<T> current;

        public:
        Progress(const std::string& name = "", ProgressEstimation estimation = ProgressEstimation::LINEAR, 
            const T initial = T(), const T objective = T());
        Progress(Progress<T>& progress);
        Progress<T>& operator=(const Progress<T>& progress);
        std::string report(const unsigned int level) const override;

        Progress<T>& operator=(const T& val);
        Progress<T>& operator+=(const T& val);
        Progress<T>& operator++();
        Progress<T>& operator-=(const T& val);
        Progress<T>& operator--();
        Progress<T>& operator*=(const T& val);
        Progress<T>& operator/=(const T& val);
        Progress<T>& operator%=(const T& val);

        operator std::atomic<T>&();
        operator T() const;
        
    };
}

#endif // PROGRESS_HPP
