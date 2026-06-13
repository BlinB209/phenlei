#ifdef WIN32
inline double TimeSpan::GetWallTime()
{
    return static_cast<double>(clock());
}

inline double TimeSpan::GetSecondTime(const double &time_elapsed)
{ 
    return time_elapsed / CLOCKS_PER_SEC;
}

inline void TimeSpan::ResetTime()
{
    last_time = GetWallTime();
    now_time = last_time;
}

inline double TimeSpan::GetTimeSpan()
{
    now_time = GetWallTime();
    double time_elapsed = now_time - last_time;
    last_time = now_time;
    return GetSecondTime(time_elapsed);
}

inline void TimeSpan::ShowTimeSpan(const string &title)
{
    now_time = GetWallTime();
    double time_elapsed = now_time - last_time;
    last_time = now_time;
    cout.setf(ios::fixed);
    cout.precision(2);
    cout << title << "time elapsed : " << GetSecondTime(time_elapsed) << " seconds" << "\n";
}

inline void TimeSpan::ShowTimeSpanToLogFile(const string &title)
{
    ostringstream oss;
    now_time = GetWallTime();
    double time_elapsed = now_time - last_time;
    last_time = now_time;

    //oss << title << " time elapsed : " << GetSecondTime(time_elapsed) << " seconds" << "\n";
    oss << setiosflags(ios::fixed) << setprecision(3) << setw(6) << GetSecondTime(time_elapsed) << "  seconds: " << title;
    WriteLogFile(oss.str());
}
#else
inline void TimeSpan::ResetTime()
{
    gettimeofday(&last_time, NULL);
}

inline double TimeSpan::GetTimeSpan()
{
    gettimeofday(&now_time, NULL);
    double time_elapsed = (now_time.tv_sec - last_time.tv_sec) + (now_time.tv_usec - last_time.tv_usec) * 1.0e-6;
    double timespan = ((now_time.tv_sec - last_time.tv_sec) + (now_time.tv_usec - last_time.tv_usec) * 1.0e-6);
    last_time = now_time;
    return timespan;
}

inline void TimeSpan::ShowTimeSpan(const string &title)
{
    cout.setf(ios::fixed);
    cout.precision(2);
    cout << title << "time elapsed : " << GetTimeSpan() << " seconds" << "\n";
}

inline void TimeSpan::ShowTimeSpanToLogFile(const string &title)
{
    PH_BarrierSepMode();

    ostringstream oss;
    oss << title << " time elapsed : " << GetTimeSpan() << " seconds" << "\n";
    WriteLogFile(oss.str());
}
#endif