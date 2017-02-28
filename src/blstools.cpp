/***************************************************************************
 *   Copyright (C) 2017 Jan Fostier (jan.fostier@ugent.be)                 *
 *   This file is part of BLStools                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <string>
#include <iostream>

#include "blstools.h"
#include "bls.h"
#include "pwmscan.h"

using namespace std;

void printProgramVersion()
{
        cout << "BLStools version " << BLSTOOLS_MAJOR_VERSION << "."
             << BLSTOOLS_MINOR_VERSION << "." << BLSTOOLS_PATCH_LEVEL << "\n";

        cout << "Copyright (C) 2017 Jan Fostier (jan.fostier@ugent.be)\n";
        cout << "This is free software; see the source for copying conditions. "
                "There is NO\nwarranty; not even for MERCHANTABILITY or "
                "FITNESS FOR A PARTICULAR PURPOSE.\n" << endl;
}

void printUsage()
{
        cout << "Usage: blstools command [options]\n\n";

        cout << " command\n";
        cout << "  scan\t\t\tscan for pwm occurrences\n";
        cout << "  bls\t\t\tcompute the branch length score\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\tdisplay help page\n";
        cout << "  -v\t--version\tdisplay version\n\n";

        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n";
}

void runBLSModule(int argc, char **argv)
{
        try {
                BLS bls(argc, argv);
        } catch (runtime_error e) {
                cerr << e.what() << endl;
                exit(EXIT_FAILURE);
        }
}

void runScanModule(int argc, char **argv)
{
        try {
                PWMScan scan(argc, argv);
        } catch (runtime_error e) {
                cerr << e.what() << endl;
                exit(EXIT_FAILURE);
        }
}

int main(int argc, char **argv)
{
        Command command = Command::none;

        // parse first parameter
        if (argc > 1) {
                string arg(argv[1]);

                if ((arg == "-h") || (arg == "--help")) {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if ((arg == "-v") || (arg == "--version")) {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                } else if (arg == "bls") {
                        command = Command::bls;
                } else if (arg == "scan") {
                        command = Command::scan;
                }
        }

        switch (command)
        {
                case Command::none:
                        printUsage();
                        exit(EXIT_FAILURE);
                        break;
                case Command::bls:
                        runBLSModule(argc, argv);
                        break;
                case Command::scan:
                        runScanModule(argc, argv);
                        break;
        }

        return EXIT_SUCCESS;
}
