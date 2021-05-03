#pragma once

#include <chrono>

struct KTimer
{
   KTimer(const char *sName, float fThresholdMS = 0) : m_sName(sName), m_fThresholdMS(fThresholdMS)
   {
       m_start = std::chrono::high_resolution_clock::now();
   }
   ~KTimer()
   {
       auto end = std::chrono::high_resolution_clock::now();
       double durationMS = std::chrono::duration_cast<std::chrono::duration<double>>(end - m_start).count() * 1000;
       if (durationMS > m_fThresholdMS)
       {
           FILE* fp = nullptr;
           fopen_s(&fp, "C:\\GitHub\\myBiology\\neuronApp\\timerOutput.txt", "a+");
           if (fp)
           {
               fprintf(fp, "%s: %.2lf ms\n", m_sName, durationMS);
               fclose(fp);
           }
       }
   }
private:
    const char* m_sName;
    float m_fThresholdMS;
    std::chrono::high_resolution_clock::time_point m_start;
};