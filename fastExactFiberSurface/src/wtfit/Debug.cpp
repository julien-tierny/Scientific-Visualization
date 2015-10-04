/*
 * file:                Debug.cpp
 * description:         Minimalist debug handling.
 * author:              (©) Julien Tierny <julien.tierny@lip6.fr>.
 * date:                February 2011.
 */

#include                <Debug.h>

bool wtfit::welcomeMsg_ = false;
bool wtfit::goodbyeMsg_ = true;
int wtfit::globalDebugLevel_ = 0;

Debug::Debug() { 
  
  debugLevel_ = 3; 
  threadNumber_ = 1;
  lastObject_ = false;
#ifdef withOpenMP
  threadNumber_ = omp_get_num_procs();
#endif
  wrapper_ = NULL;
  
  // avoid warnings
  if(goodbyeMsg_) goodbyeMsg_ = true; 
  
  if(welcomeMsg_){
    welcomeMsg_ = false;
    stringstream s;
    s << "[Common] "
"__        _______ _____ ___ _____        __  __    ____   ___  _ ____"
            << endl << "[Common] "
"\\ \\      / /_   _|  ___|_ _|_   _|      / /__\\ \\  |___ \\ / _ \\/ | ___|"
            << endl << "[Common] "
" \\ \\ /\\ / /  | | | |_   | |  | |       | |/ __| |   __) | | | | |___ \\"
            << endl << "[Common] "
"  \\ V  V /   | | |  _|  | |  | |       | | (__| |  / __/| |_| | |___) |"
            << endl << "[Common] "
"   \\_/\\_/    |_| |_|   |___| |_|       | |\\___| | |_____|\\___/|_|____/"
            << endl << "[Common] "
"                                        \\_\\  /_/" << endl;
    s << "[Common] Welcome!" << endl;
    dMsg(cout, s.str(), 1);
  }
}

Debug::~Debug(){
  if((lastObject_)&&(goodbyeMsg_)){
    stringstream msg;
    msg << "[Common] Goodbye :)" << endl;
    dMsg(cout, msg.str(), 1);
    goodbyeMsg_ = false;
  }
}

const int Debug::dMsg(ostream &stream, string msg, 
  const int &debugLevel) const{
  if((debugLevel_ >= debugLevel)||(globalDebugLevel_ >= debugLevel))
    stream << msg.data();
  return 0;
}

const int Debug::err(const string msg, const int &debugLevel) const{
  return dMsg(cerr, msg, 0);
}

const int Debug::msg(const char *msg, const int &debugLevel) const{
  return dMsg(cout, string(msg), debugLevel);
}

const int Debug::setDebugLevel(const int &debugLevel){
  debugLevel_ = debugLevel;
  return 0;
}

int Debug::setWrapper(Wrapper *wrapper, const int &debugLevel){
  wrapper_ = wrapper;
  debugLevel_ = debugLevel;
  return 0;
}