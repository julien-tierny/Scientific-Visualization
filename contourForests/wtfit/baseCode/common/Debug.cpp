#include                <Debug.h>

bool wtfit::welcomeMsg_ = true;
bool wtfit::goodbyeMsg_ = true;
int wtfit::globalDebugLevel_ = -INT_MAX;

Debug::Debug() { 
  
  threadNumber_ = 1;
  lastObject_ = false;
#ifdef withOpenMP
  threadNumber_ = omp_get_num_procs();
#endif
  wrapper_ = NULL;
 
  debugLevel_ = infoMsg;
  
  // avoid warnings
  if(goodbyeMsg_) goodbyeMsg_ = true; 
  
  if((wtfit::welcomeMsg_)&&(debugLevel_)){
    wtfit::welcomeMsg_ = false;
    stringstream s;
      
    s << "[Common] "
"__        _______ _____ ___ _____        __  __    ____   ___  _  __"
      << endl << "[Common] "
"\\ \\      / /_   _|  ___|_ _|_   _|      / /__\\ \\  |___ \\ / _ \\/ |/ /_"
      << endl << "[Common] "
" \\ \\ /\\ / /  | | | |_   | |  | |       | |/ __| |   __) | | | | | '_ \\"
      << endl << "[Common] "
"  \\ V  V /   | | |  _|  | |  | |       | | (__| |  / __/| |_| | | (_) |"
      << endl << "[Common] "
"   \\_/\\_/    |_| |_|   |___| |_|       | |\\___| | |_____|\\___/|_|\\___/"
      << endl << "[Common] "
"                                        \\_\\  /_/"
      << endl;
    s << "[Common] Welcome!" << endl;
    dMsg(cout, s.str(), 1);
  }
}

Debug::~Debug(){
  if((lastObject_)&&(wtfit::goodbyeMsg_)){
    stringstream msg;
    msg << "[Common] Goodbye :)" << endl;
    dMsg(cout, msg.str(), 1);
    wtfit::goodbyeMsg_ = false;
  }
}

const int Debug::dMsg(ostream &stream, string msg, 
  const int &debugLevel) const{
 
  if(!debugLevel_) 
    return -1;
  
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

int Debug::setWrapper(const Wrapper *wrapper){
  
  wrapper_ = (Wrapper *) wrapper;
  setDebugLevel(((Debug *) wrapper)->debugLevel_);
  setThreadNumber(((Debug *) wrapper)->threadNumber_);
  return 0;
}
