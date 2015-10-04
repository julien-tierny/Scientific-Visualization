/*
 * file:                Debug.h
 * description:         Minimalist debug handling.
 * author:              (©) Julien Tierny <tierny@telecom-paristech.fr>.
 * date:                February 2011.
 */

#ifndef                 _DEBUG_H
#define                 _DEBUG_H

#include                <fstream>
#include                <iostream>
#include                <sstream>
#include                <string>
#include                <vector>

#ifdef __linux__
#include                <unistd.h>
#include                <sys/time.h>
#endif

using namespace std;

class Debug{

  public:

    Debug() { debugLevel_ = 2; };
    ~Debug(){};

    inline const int dMsg(ostream &stream, string msg, 
      const int &debugLevel) const{
      if(debugLevel_ >= debugLevel)
        stream << msg.data();
      return 0;
    }

    virtual inline const int setDebugLevel(const int &debugLevel){
      debugLevel_ = debugLevel;
    };

  protected:

    int                     debugLevel_;
};

class DebugTimer : public Debug{

  public:

    DebugTimer(){ 
#ifdef __linux__
      gettimeofday(&start_, NULL);
#endif
    };

    DebugTimer(const DebugTimer &other){ 
#ifdef __linux__
      start_ = other.start_;
#endif
    }

    inline double getElapsedTime(){
      
#ifdef __linux__
      struct timeval end;
      gettimeofday(&end, NULL);

      return ((end.tv_sec*1000000 + end.tv_usec)
          - (start_.tv_sec*1000000 
            + start_.tv_usec))/1000000.0;
#else
      return 0;
#endif
    };

    inline double getTime(){ 

#ifdef __linux__
      return (start_.tv_sec*1000000 + start_.tv_usec)/1000000.0;
#else
      return 0;
#endif
      
    };

    inline void reStart(){ 
#ifdef __linux__
      gettimeofday(&start_, NULL);
#endif
    }

  protected:

#ifdef __linux__
    struct timeval  start_;
#endif
};

class DebugMemory : public Debug{
  
  public:
    
    inline float getUsage(){
#ifdef __linux__
      // horrible hack since getrusage() doesn't seem to work well under linux
      stringstream procFileName;
      procFileName << "/proc/" << getpid() << "/statm";
      
      ifstream procFile(procFileName.str().data(), ios::in);
      if(procFile){
        float memoryUsage;
        procFile >> memoryUsage;
        procFile.close();
        return memoryUsage/1024.0;
      }
      
      return 0;
#endif
    };
    
  protected:
};

#endif
