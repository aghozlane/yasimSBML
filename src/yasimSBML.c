/**
 * This file is part of yasimSBML (http://www.labri.fr/perso/ghozlane/metaboflux/about/yasimSBML.php)
 * Copyright (C) 2010 Amine Ghozlane from LaBRI and University of Bordeaux 1
 *
 * yasimSBML is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * yasimSBML is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file yasimSBML.c
 * \brief Main program
 * \author {Amine Ghozlane}
 * \version 1.0
 * \date 27 octobre 2009
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <sbml/SBMLTypes.h>
#include "especes.h"
#include "simulation.h"

/**
 * \fn int main(int argc, char **argv)
 * \author Amine Ghozlane
 * \brief  Enter in the program for simulated annealing
 * \param  argc Number of arguments
 * \param  argv List of arguments
 * \return EXIT_SUCCESS Normal stop of the program
 */
int main(int argc, char **argv) {
  SBMLDocument_t *file=NULL;
  Model_t *mod=NULL;
  unsigned int errors = 0;
  /*unsigned int level;
  unsigned int version;*/

  /* File adress */
  if (argc <= 2) {
      fprintf(stderr, "Usage ./yasimSBML.exe sbml_file.xml number_of_simulation\n");
      return 1;
  }

  /* Load File  */
  file = readSBML(argv[1]);
  errors = SBMLDocument_getNumErrors(file);

  /* Error File */
  if (errors > 0) {
      printf("Encountered the following SBML error(s):\n");
      SBMLDocument_printErrors(file, stdout);
      printf("Simulation skipped.  Please correct the problems above first.\n");
      return errors;
  }

  /*level = SBMLDocument_getLevel(file);
  version = SBMLDocument_getVersion(file);
  printf("SBML Level : %d, version : %d\n", level, version);*/

  /* Testing mod file */
  mod = SBMLDocument_getModel(file);
  if (mod == NULL)
    fprintf(stderr, "Error: file %s doesn't exist!\n", argv[1]);
  /*compute_simulation(mod);*/
  SBML_compute_simulation_mean(mod,atoi(argv[2]));

  /* free memory */
  Model_free(mod);
  /*SBMLDocument_free(file);*/

  return EXIT_SUCCESS;
}

