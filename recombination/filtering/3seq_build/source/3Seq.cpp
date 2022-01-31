/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   3Seq.cpp
 * CREATED ON:  02 June 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: The main source file for executable 3SEQ.
 *              This file contains main() function.
 *
 * HISTORY:     Version     Date            Description
 *              1.1104      2011-04-07      Original source from 3SEQ version 1.110407
 *              1.1106      2011-06-02      Reorganise the whole program (follow OOP design).
 *
 * VERSION:     1.1106
 * LAST EDIT:   02 June 2011
 */



//#define TESTING /* Use this flag whenever you want to run unit tests */
#ifdef TESTING

#include <iostream>

#include "AppInfo.h"
#include "run-mode/Run.h"
#include "run-mode/FullRun.h"

int main(int argc, char ** argv) {
    int _argc = 2;
    std::string argList[] = {"3seq", "-f", "vvv", "include.mtDNA.aln"};
    
    char ** _argv = new char * [_argc];
    
    for (int i = 0; i < _argc; i++) {
        _argv[i] = const_cast<char *> (argList[i].c_str());
    }
    
    if (_argc >= 1) {
        /* Store the path to the executable file. */
        AppInfo::instance()->setExecutablePath(_argv[0]);
    }
    std::cout << AppInfo::instance()->getVersionFullString() << std::endl;
    
    srand(time(nullptr));
    for (int i = 0; i < 10; i++) {
        std::cout << rand() << std::endl;
    }
    
    
    return 0;
}

#else


#include <iostream>
#include <cmath>
#include <ctime>

#include "AppInfo.h"
#include "Interface.h"
#include "run-mode/Run.h"

int main(int argc, char** argv) {
    /* Seed the random number generator for later uses. */
    srand(time(nullptr));
    
    if (argc >= 1) {
        /* Store the path to the executable file. */
        AppInfo::instance()->setExecutablePath(argv[0]);
    }
    
    Run* run = nullptr;

    try {
        /* Detect run-mode and return a run-object */
        run = Run::getRun(argc, argv);

        if (run == nullptr) {
            /* Cannot detect run-mode, just show the program usage */
            Interface::instance().startProgram(AppInfo::instance()->getDescription());
            Interface::instance().showPTableLink();
            Interface::instance().catchUsageModes();
            Interface::instance().showLog(true);

        } else {
            /* Here is the main flow of every analysis (of all kinds) */
            run->parseCmdLine();
            run->initLogFile();
	//TODO: edited
            run->perform();
            //run->perform_output();
        }

        /* Throw a SUCCESS signal to indicate that the analysis has been done successfully.
         * Besides, this method is also necessary to save the log file. */
        if (run == nullptr || run->getMode() != Run::Mode::WRITE) {
            /* The write mode will throw Exit signal by itself */
            Interface::instance().throwExitSignal(0);
        }


    } catch (Interface::ExitSignal) {
        /* Do nothing */


    } catch (std::bad_alloc) {
        /* Reaching here means that no ExitSignal has been thrown and
         * the program is terminated due to out-of-memory. */
        Interface::instance() << "Out of memory.\n";
        try {
            Interface::instance().showError(true, true);
        } catch (...) {
            // Do nothing
        }


    } catch (...) {
        /* Reaching here means that no ExitSignal has been thrown and
         * the program is terminated by an unknown error. */
        Interface::instance()
                << "An unexplained exception occured during the analysis.\n"
                << "For more help, email author: " << AppInfo::instance()->getAuthorContact() << "\n";
        try {
            Interface::instance().showError(true, true);
        } catch (...) {
            // Do nothing
        }
    }

    if (run != nullptr) {
        delete run;
        run = nullptr;
    }

    return 0;
}

#endif /* TESTING */
