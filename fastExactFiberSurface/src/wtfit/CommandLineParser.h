/*
 * file:                CommandLineParser.h
 * description:         Basic command line parsing.
 * author:              (©) Julien Tierny <julien.tierny@lip6.fr>
 * date:                September 2013.
 */

#ifndef                 _COMMAND_LINE_PARSER_H
#define                 _COMMAND_LINE_PARSER_H

#include                <Debug.h>

namespace wtfit{

  class CommandLineParser : public Debug{
    
    public:
      
      // 1) constructors, destructors, operators, etc.
      class CommandLineArgument : public Debug{
        
        
        public:
          
          CommandLineArgument(){
            boolValue_ = NULL;
            intValue_ = NULL;
            realValue_ = NULL;
            stringValue_ = NULL;
          };
                  
          int print(stringstream &s) const{
            s << "[CommandLine]   ";
            if((isAnOption_)||(isOptional_)){
              s << "[";
            }  
            
            s << "-" << key_;
            
            if(isAnOption_)
              s << ":";
            
            s << " ";
            if(!isAnOption_)
              s << "<";
            
            if(description_.length()){
              s << description_;
            }
            else{
              s << "no description";
            }
            
            if(!isAnOption_)
              s << ">";
            
            if((isAnOption_)||(isOptional_)){
              s << "]";
            }
            s << endl;
            
            return 0;
          };
         
          bool              isOptional_, isAnOption_;
          bool              *boolValue_;
          int               *intValue_;
          real              *realValue_;
          string            *stringValue_;
          
          string            key_;
          string            description_;
      };
      
      CommandLineParser(){
        setIntArgument("d", &(wtfit::globalDebugLevel_),
          "Global debug level", true);
      };
      
      ~CommandLineParser(){};
      
      // 2) functions
      inline int parse(int argc, char **argv) const{
       
        for(int i = 0; i < argc; i++){
          
          if((string(argv[i]) == "-h")||(string(argv[i]) == "--help")){
            printUsage(argv[0]);
          }
          
          for(int j = 0; j < (int) arguments_.size(); j++){
            
            if(!arguments_[j].isAnOption_){
              if((string(argv[i]) == "-" + arguments_[j].key_)
                &&(i + 1 < argc)){
              
                if(arguments_[j].stringValue_){
                  // let's process a string
                  (*arguments_[j].stringValue_) = argv[i + 1];
                }
                else if(arguments_[j].intValue_){
                  stringstream s(argv[i + 1]);
                  s >> *(arguments_[j].intValue_);
                }
                else if(arguments_[j].realValue_){
                  stringstream s(argv[i + 1]);
                  s >> *(arguments_[j].realValue_);
                }
              }
            }
            else{
              if(string(argv[i]) == "-" + arguments_[j].key_){
                *(arguments_[j].boolValue_) = true;
              }
            }
          }
        }
        
        // check all the necessary arguments have been provided
        for(int i = 0; i < (int) arguments_.size(); i++){
          if(!arguments_[i].isOptional_){
            if(arguments_[i].stringValue_){
              if(!arguments_[i].stringValue_->length()){
                stringstream msg;
                msg << "[CommandLine] Missing mandatory argument:" << endl;
                arguments_[i].print(msg);
                dMsg(cerr, msg.str(), 1);
                printUsage(argv[0]);
              }
            }
            else if(arguments_[i].intValue_){
              if(*(arguments_[i].intValue_) == -INT_MAX){
                stringstream msg;
                msg << "[CommandLine] Missing mandatory argument:" << endl;
                arguments_[i].print(msg);
                dMsg(cerr, msg.str(), 1);
                printUsage(argv[0]);
              }
            }
            else if(arguments_[i].realValue_){
              if(*(arguments_[i].realValue_) == -REAL_MAX){
                stringstream msg;
                msg << "[CommandLine] Missing mandatory argument:" << endl;
                arguments_[i].print(msg);
                dMsg(cerr, msg.str(), 1);
                printUsage(argv[0]);
              }
            }
          }
        }
        
        return 0;
      };
      
      inline int printArgs(ostream &o) const{
        
        o << "[CommandLine] Options and arguments:" << endl;
        for(int i = 0; i < (int) arguments_.size(); i++){
          o << "[CommandLine]   -" << arguments_[i].key_;
          o << ": ";
          
          if(arguments_[i].isAnOption_){
            if(arguments_[i].boolValue_){
              if(*(arguments_[i].boolValue_))
                o << "true";
              else
                o << "false";
            }
            else{
              o << "(not set)";
            }
          }
          else if(arguments_[i].stringValue_){
            if(arguments_[i].stringValue_->length()){
              o << *(arguments_[i].stringValue_);
            }
            else{
              o << "(not set)";
            }
          }
          else if(arguments_[i].intValue_){
            if(*(arguments_[i].intValue_) == -INT_MAX){
              o << "(not set)";
            }
            else{
              o << *(arguments_[i].intValue_);
            }
          }
          else if(arguments_[i].realValue_){
            if(*(arguments_[i].realValue_) == -REAL_MAX){
              o << "(not set)";
            }
            else{
              o << *(arguments_[i].realValue_);
            }
          }
          
          o << endl;
        }
        
        return 0;
      }
      
      inline int printUsage(const string &binPath) const {
        
        stringstream msg;
        msg << "[CommandLine]" << endl;
        msg << "[CommandLine] Usage:" << endl;
        msg << "[CommandLine]   " << binPath << endl;
        msg << "[CommandLine] Argument(s):" << endl;
        for(int i = 0; i < (int) arguments_.size(); i++){
          if(!arguments_[i].isAnOption_){
            arguments_[i].print(msg);
          }
        }
        msg << "[CommandLine] Option(s):" << endl;
        for(int i = 0; i < (int) arguments_.size(); i++){
          if(arguments_[i].isAnOption_){
            arguments_[i].print(msg);
          }
        }
        
        dMsg(cerr, msg.str(), 1);
        
        exit(0);
        return 0;
      };
      
      inline int setOption(const string &key, bool *value,
        const string &description = ""){
        
        arguments_.resize(arguments_.size() + 1);
        arguments_.back().isOptional_ = true;
        arguments_.back().key_ = key;
        arguments_.back().description_ = description;
        arguments_.back().boolValue_ = value;
        *(arguments_.back().boolValue_) = false;
        arguments_.back().isAnOption_ = true;
       
        return 0;
      };
      
      inline int setRealArgument(const string &key, real *value, 
        const string &description = "",
        const bool &optional = false){
        
        arguments_.resize(arguments_.size() + 1);
        arguments_.back().isOptional_ = optional;
        arguments_.back().key_ = key;
        arguments_.back().description_ = description;
        arguments_.back().realValue_ = value;
        arguments_.back().isAnOption_ = false;
       
        *(arguments_.back().realValue_) = -REAL_MAX;
        
        return 0;
      };
      
      inline int setIntArgument(const string &key, int *value, 
        const string &description = "",
        const bool &optional = false){
        
        arguments_.resize(arguments_.size() + 1);
        arguments_.back().isOptional_ = optional;
        arguments_.back().key_ = key;
        arguments_.back().description_ = description;
        arguments_.back().intValue_ = value;
        arguments_.back().isAnOption_ = false;
       
        *(arguments_.back().intValue_) = -INT_MAX;
        
        return 0;
      };
      
      int setStringArgument(const string &key, string *value, 
        const string &description = "",
        const bool &optional = false){
     
        arguments_.resize(arguments_.size() + 1);
        arguments_.back().isOptional_ = optional;
        arguments_.back().key_ = key;
        arguments_.back().description_ = description;
        arguments_.back().stringValue_ = value;
        arguments_.back().isAnOption_ = false;
       
        *(arguments_.back().stringValue_) = "";
        
        return 0;
      };
      
      
    protected:
      
      vector<CommandLineArgument>
                        arguments_;
  };
}

#endif
