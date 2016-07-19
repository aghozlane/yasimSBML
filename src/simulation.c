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
 * \file simulation.c
 * \brief Simulate a petri net
 * \author {Amine Ghozlane}
 * \version 1.0
 * \date 27 octobre 2009
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <sbml/SBMLTypes.h>
#include "especes.h"
#include "simulation.h"

/*TODO Verification des BoundaryConditions,Local parameter values, Global parameter values, Reaction, des equations */

/**
 * \fn void SBML_initEspeceAmounts(Model_t *mod, pEspeces molecules, int nbEspeces)
 * \author Amine Ghozlane
 * \brief  Alloc memory and initialize the struct Especes
 * \param  mod Model of the SBML file
 * \param  molecules Struct Especes
 * \param  nbEspeces Number of molecules
 */
void SBML_initEspeceAmounts(Model_t *mod, pEspeces molecules, int nbEspeces) {
  int i;
  Species_t *esp;

  /* Initialisation des quantites des especes*/
  for (i = 0; i < nbEspeces; i++) {
      esp = Model_getSpecies(mod, i);
      /*printf("espece: %s, compartiment : %s\n",Species_getId(esp),Species_getCompartment(esp));*/
      Especes_save(molecules, i, Species_getInitialAmount(esp),
          Species_getId(esp), Species_getCompartment(esp));
  }
}

/**
 * \fn void SBML_setReactions(Model_t *mod, pEspeces molecules, pScore result, double *reactions_ratio, int nbReactions, int nbEspeces)
 * \author Amine Ghozlane
 * \brief  Alloc memory and initialize the struct Especes
 * \param  mod Model of the SBML file
 * \param  molecules Struct Especes
 * \param  result Struct Score
 * \param  reactions_ratio List of computed reaction ratio
 * \param  nbReactions Number of reaction
 * \param  nbEspeces Number of molecules
 */
void SBML_setReactions(Model_t *mod, pEspeces molecules, int nbReactions, int nbEspeces) {
  int ref = 0, i, j;
  SpeciesReference_t *reactif;
  Species_t *especeId;
  Reaction_t *react;
  const char *kf;
  const ASTNode_t *km;
  KineticLaw_t *kl;

  /*fprintf(stdout, "Debut reaction\n");*/
  /* Recherche les reactions ou apparaissent chaque espece */
  for (i = 0; i < nbReactions; i++) {
      react = Model_getReaction(mod, i);
      for (j = 0; j < (int)Reaction_getNumReactants(react); j++) {
          reactif = Reaction_getReactant(react, j);
          especeId = Model_getSpeciesById(mod, SpeciesReference_getSpecies(
              reactif));
          ref = Especes_find(molecules, Species_getId(especeId), nbEspeces);
          /*printf("Reactif :Quantite de ref %d :%f \n", ref,Especes_getQuantite(molecules, ref));*/

          kl = Reaction_getKineticLaw(react);

          if (KineticLaw_isSetFormula(kl)) {
              kf = KineticLaw_getFormula(kl);
              /*printf("c du kf %s\n", kf);*/
              Especes_allocReactions(molecules, ref, react,
                  SBML_evalExpression(kf));
          } else {
              km = KineticLaw_getMath(kl);
              printf("C du kM... cas encore ignore\n");
              exit(1);
          }
      }
  }
  /*fprintf(stdout, "Fin reaction\n");*/
}

/**
 * \fn double SBML_evalExpression(const char *formule)
 * \author Amine Ghozlane
 * \brief  Get the reaction ratio define in the sbml
 * \param  formule Formule SBML
 * \return Return double value of the constraint
 */
double SBML_evalExpression(const char *formule) {
  return atof(formule);
}

/**
 * \fn int SBML_checkQuantite(Model_t *mod, Reaction_t *react, int nbEspeces, pEspeces molecules)
 * \author Amine Ghozlane
 * \brief  Determine the number of reaction for one molecule
 * \param  mod Model of the SBML file
 * \param  react Reaction id
 * \param  nbEspeces Number of molecules
 * \param  molecules Struct Especes
 * \return Number of reaction for one molecule
 */
int SBML_checkQuantite(Model_t *mod, Reaction_t *react, int nbEspeces, pEspeces molecules)
{
  double quantite = 0.0, minStep = 0.0, temp = 0.0;
  int ref = 0, i;
  SpeciesReference_t *reactif;
  Species_t *especeId;

  reactif = Reaction_getReactant(react, 0);
  especeId = Model_getSpeciesById(mod, SpeciesReference_getSpecies(reactif));
  ref = Especes_find(molecules, Species_getId(especeId), nbEspeces);

  if ((quantite = Especes_getQuantite(molecules, ref)) <= 0.0) {
      /*printf("je fais le return\n");*/
      return END;
  }
  /*printf("quantite : %d\n",quantite);*/
  minStep = quantite / SpeciesReference_getStoichiometry(reactif);
  /*printf("quantite : %d, minStep init : %d\n",quantite,minStep);*/
  /* Cas ou le nombre de pas est egal a 0 */
  if(minStep==0.0) return END;

  for (i = 1; i < (int)Reaction_getNumReactants(react); i++) {
      reactif = Reaction_getReactant(react, i);
      especeId = Model_getSpeciesById(mod, SpeciesReference_getSpecies(reactif));
      ref = Especes_find(molecules, Species_getId(especeId), nbEspeces);
      quantite= Especes_getQuantite(molecules, ref);
      /* La quantite est egale a 0 */
      if(quantite<=0.0) return END;
      temp = floor(Especes_getQuantite(molecules, ref)/SpeciesReference_getStoichiometry(reactif));
      /* Cas ou le nombre de pas est egal a 0 */
      if(temp==0.0) return END;
      /* Si le nouveau nombre est inferieur au precedent, on change la valeur de minStep */
      if (minStep > temp) minStep = temp;
  }

  return (int)minStep;
}

/**
 * \fn Reaction_t * SBML_reactChoice(pEspeces molecules, const gsl_rng * r, int ref)
 * \author Amine Ghozlane
 * \brief  Determine randomly the reaction to achieve for several nodes reactions
 * \param  molecules Struct Especes
 * \param  r Random number generator
 * \param  ref Number reference of one molecule
 * \return Id of the selected reaction
 */
Reaction_t * SBML_reactChoice(pEspeces molecules, const gsl_rng * r, int ref) {

  pReaction temp = NULL;
  pReaction Q = molecules[ref].system;
  double value = gsl_rng_uniform(r) * 100.0;
  double choice = 0.0;

  /*printf("Choix de la reaction :\n");
	 printf("value : %f\n", value);*/

  /* Choix de la reaction a realiser */
  if (value < Q->ratio)
    temp = Q;
  else {
      do {
          /*printf("reaction %s : ratio %f\n", Reaction_getId(Q->link), Q->ratio);*/
          choice += Q->ratio;
          /*printf("choice : %f\n", choice);*/
          Q = Q->suivant;
          temp = Q;
          /*printf("ratio suivant : %f\n",Q->ratio );*/
      } while (Q->suivant != NULL && value <= (choice + Q->suivant->ratio)
          && value > choice);
  }
  if (temp == NULL) {
      fprintf(stderr, "on a un probleme de ratio\n");
      exit(EXIT_FAILURE);
  }

  /*printf("Resultat :\n");
	 printf("value : %f, choice :%f\n", value, (choice + Q->ratio));
	 printf("reaction %s : ratio %f\n", Reaction_getId(temp->link), temp->ratio);*/
  /*if (Q->suivant == NULL)
	 printf("choice :%f", choice);*/
  /*printf("c bon\n");*/

  return (temp->link);

}

/**
 * \fn void SBML_reaction(Model_t *mod, pEspeces molecules, Reaction_t *react, int nbEspeces)
 * \author Amine Ghozlane
 * \brief  Simulation of a discrete transision
 * \param  mod Model of the SBML file
 * \param  molecules Struct Especes
 * \param  react Reaction id
 * \param  nbEspeces Number of molecules
 */
void SBML_reaction(Model_t *mod, pEspeces molecules, Reaction_t *react, int nbEspeces)
{
  SpeciesReference_t *reactif;
  Species_t *especeId;
  int i, ref = 0;

  /*boucle pour retirer des reactifs*/
  for (i = 0; i < (int)Reaction_getNumReactants(react); i++) {
      reactif = Reaction_getReactant(react, i);
      especeId = Model_getSpeciesById(mod, SpeciesReference_getSpecies(reactif));
      ref = Especes_find(molecules, Species_getId(especeId), nbEspeces);

      /*printf("Reactif : %s",Species_getId(especeId) );
		 printf("ref : %d\n",ref);
		 printf("Reactif :Quantite de ref %d :%d \n",ref,Espece_getQuantite(molecules,ref));*/
      Especes_setQuantite(molecules, ref,(Especes_getQuantite(molecules, ref)- SpeciesReference_getStoichiometry(reactif)));
      /*printf("Apres Reactif: Quantite de ref %d :%d \n",ref,Espece_getQuantite(molecules,ref));*/
  }

  /*boucle pour ajouter des produits */
  for (i = 0; i < (int)Reaction_getNumProducts(react); i++) {
      reactif = Reaction_getProduct(react, i);
      especeId = Model_getSpeciesById(mod, SpeciesReference_getSpecies(reactif));
      /*printf("Produit : %s",Species_getId(especeId) );*/
      ref = Especes_find(molecules, Species_getId(especeId), /*Species_getCompartment(especeId),*/ nbEspeces);
      /*printf("Produit : ref : %d\n",ref);
		 printf("Produit : Quantite de ref %d :%d \n",ref,Espece_getQuantite(molecules,ref));*/

      Especes_setQuantite(molecules, ref,(Especes_getQuantite(molecules, ref)+ SpeciesReference_getStoichiometry(reactif)));
      /*printf("Apres Produit : Quantite de ref %d :%d \n",ref,Espece_getQuantite(molecules,ref));*/
  }
}

/**
 * \fn void SBML_allocTest(pTestReaction T, int nbReactions)
 * \author Amine Ghozlane
 * \brief  Alloc memory and initialize the struct pTestReaction
 * \param  T Empty struct TestReaction
 * \param  nbReactions Number of reactions
 */
void SBML_allocTest(pTestReaction T, int nbReactions)
{
  int i;

  /* Initialisation des tableaux des reactions */
  T->tabReactions= (Reaction_t **) malloc(nbReactions * sizeof(Reaction_t *));
  T->minStepTab = (int*) malloc(nbReactions * sizeof(int));

  if (T->tabReactions == NULL)
    exit(EXIT_FAILURE);
  if (T->minStepTab == NULL)
    exit(EXIT_FAILURE);

  for (i = 0; i < nbReactions; i++) {
      T->tabReactions[i] = NULL;
      T->minStepTab[i] = 0;
  }
}

/**
 * \fn void SBML_freeTest(pTestReaction T)
 * \author Amine Ghozlane
 * \brief  Free memory of the struct TestReaction
 * \param T Struct TestReaction gives data on reaction
 */
void SBML_freeTest(pTestReaction T)
{
  if (T->tabReactions != NULL)
    free(T->tabReactions);
  if (T->minStepTab != NULL)
    free(T->minStepTab);
}

/**
 * \fn int SBML_EstimationReaction(Model_t *mod, pTestReaction T, pEspeces molecules, int ref, int nbEspeces)
 * \author Amine Ghozlane
 * \brief  Alloc memory and initialize the struct Especes
 * \param  mod Model of the SBML file
 * \param  T Struct TestReaction gives data on reaction
 * \param  molecules Struct Especes
 * \param  ref Number reference of one molecule
 * \param  nbEspeces Number of molecules
 * \return Estimated number of feasible step by reaction
 */
int SBML_EstimationReaction(Model_t *mod, pTestReaction T, pEspeces molecules, int ref, int nbEspeces)
{
  pReaction Q = molecules[ref].system;
  int i = 0, min = 0, curr = 0;

  /* Compte le nombre de reactions rattachees a une espece */
  while (Q != NULL) {
      T->tabReactions[i] = Q->link;
      T->minStepTab[i] = SBML_checkQuantite(mod, Q->link, nbEspeces,
          molecules);
      if (i == 0)
        min = T->minStepTab[i];
      curr = T->minStepTab[i];
      if (T->minStepTab[i] <= 0)
        return END;
      else if (curr < min)
        min = curr;
      Q = Q->suivant;
      i++;
  }

  return min;
}

/**
 * \fn int SBML_simulate(Model_t *mod, pEspeces molecules, const gsl_rng * r, pTestReaction T, char **banned, int nbBanned, int nbEspeces, int ref)
 * \author Amine Ghozlane
 * \brief  Simulate one step of petri net
 * \param  mod Model of the SBML file
 * \param  molecules Struct Especes
 * \param  r Random number generator
 * \param  T Struct TestReaction gives data on reaction
 * \param  banned List of banned compound
 * \param  nbBanned Number of banned compound
 * \param  nbEspeces Number of molecules
 * \param  ref Number reference of one molecule
 * \return Condition of stop/pursue
 */
int SBML_simulate(Model_t *mod, pEspeces molecules, const gsl_rng * r, pTestReaction T, int nbEspeces, int ref)
{
  Reaction_t *react = NULL;
  int minStep = 0, valid = 0/*, i=0*/;
  /*int nbReactions = Especes_getNbreactions(molecules, ref);*/
  int nbReactions = molecules[ref].nbReactions;

  /*printf("nbreactions : %d, ref : %d\n", nbReactions, ref);
	 printf("molecules : %s\n", molecules[ref].id);*/

  /* Probleme ATP, ADP, NADH, NAD+ */
  if (!strcmp(molecules[ref].id, "ADP") || !strcmp(molecules[ref].id, "ATP")
      || !strcmp(molecules[ref].id, "NADH") || !strcmp(molecules[ref].id,
          "NADplus") || !strcmp(molecules[ref].id, "NADPH") || !strcmp(
              molecules[ref].id, "NADPplus")|| !strcmp(molecules[ref].id, "NADplus_g")||
              !strcmp(molecules[ref].id, "NADH_g")|| !strcmp(molecules[ref].id, "ADP_g")||
              !strcmp(molecules[ref].id, "ATP_g")|| !strcmp(molecules[ref].id, "ADP_c")|| !strcmp(molecules[ref].id, "ATP_c")) {
      /*printf("END : molecules : %s\n", molecules[ref].id);*/
      return END;
  }

  /* Elimine les calculs sur ADP... */
  /*for(j=0;j<BANNI;j++){
	 if(!strcmp(tab[i],molecules[ref].id)){
	 printf("test : tab: %s molecules : %s\n",tab[i],molecules[ref].id );
	 return END;
	 }
	 }*/

  /* Variation selon le cas. */
  switch (nbReactions) {
  /* Cas ou aucune reaction n'est possible */
  case 0:
    /*printf("Cas N°1\n");*/
    return END;
    break;
    /* Cas ou une seule reaction est possible */
  case 1:
    /*printf("Cas N°2\n");*/
    react = molecules[ref].system->link;
    if ((minStep = SBML_checkQuantite(mod, react, nbEspeces, molecules))<= END)
      return END;
    /*printf("minstep = %d\n", minStep);*/

    while (minStep > 0) {
        SBML_reaction(mod, molecules, react, nbEspeces);
        minStep--;
    }
    break;
    /* Cas ou  plusieurs reactions sont possibles*/
  default:
    /*printf("Cas N°3\n");*/
    /* Allocation de la memoire au tableau des reactions */
    SBML_allocTest(T, nbReactions);
    valid = SBML_EstimationReaction(mod, T, molecules, ref, nbEspeces);
    if (valid <= END) {
        /*printf("on ne fait rien\n");*/
        return END;
    }
    while (Especes_getQuantite(molecules, ref) > 0.0 && valid > END) {
        react = SBML_reactChoice(molecules, r, ref);
        /*printf("molecule : %s, reaction %s\n", molecules[ref].id, Reaction_getId(react));*/
        /*printf("reaction : %s\n", Reaction_getId(react));*/
        /*printf("valid : %d, i : %d\n", valid, i);*/
        SBML_reaction(mod, molecules, react, nbEspeces);
        /*if(i>=valid) valid=SBML_EstimationReaction(mod, T, molecules, ref, nbEspeces);*/
        /*i++;*/
        valid--;
    }
    /*printf("Fin de simulation\n");*/

    /* Liberation de la memoire allouee au tableau des reactions */
    SBML_freeTest(T);
    break;
  }
  /*printf("On quitte la simulation\n");*/
  return PURSUE;
}

/**
 * \fn void SBML_compute_simulation(pScore result, Model_t *mod, double *reactions_ratio, gsl_rng * r, char **banned, int nbBanned)
 * \author Amine Ghozlane
 * \brief  Simulation of metabolic network
 * \param  result Struct Score
 * \param  mod Model of the SBML file
 * \param  reactions_ratio List of computed reaction ratio
 * \param  r Random number generator
 * \param  banned List of banned compound
 * \param  nbBanned Number of banned compound
 */
void compute_simulation(Model_t *mod)
{
  int i, nbReactions = 0, nbEspeces = 0, temp = 1, tempo = 0;
  pEspeces molecules;
  const gsl_rng_type *T;
  gsl_rng * r;
  pTestReaction TR = NULL;

  /* Verifications */
  /*TODO numGlobalParameters=Model_getNumParameters(mod);
	checkBoundaryConditions(mod);
	if(numGlobalParameters>0){
	    globalParVec=gsl_vector_alloc(numGlobalParameters);
	    setGlobalParVec();
	  }*/
  /* Allocation memoire et initialisation des generateurs de nombre aleatoire */
  gsl_rng_env_setup();
  if (!getenv("GSL_RNG_SEED"))
    gsl_rng_default_seed = time(NULL);
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_default_seed = gsl_rng_uniform(r);

  /* Allocation memoire */
  TR = (pTestReaction) malloc(1 * sizeof(TestReaction));
  if (T == NULL)
    exit(EXIT_FAILURE);

  /* File information */
  nbReactions = Model_getNumReactions(mod);
  nbEspeces = Model_getNumSpecies(mod);
  /*printf("Nom du model : %s\n",Model_getId(mod));*/
  molecules = Especes_alloc(nbEspeces);
  /* Initialisation de la quantite des especes */
  SBML_initEspeceAmounts(mod, molecules, nbEspeces);
  /* Initialisation des reactions et des ratios*/
  SBML_setReactions(mod, molecules, nbReactions, nbEspeces);


  /*TODO checkRatio(mod,molecules);*/

  /* Test des donnees enregistrees */
  printf("\nEtat des especes :\n\n");
  Especes_print(molecules, nbEspeces);

  printf("\nDebut de simulation ...\n");
  /* Simulation des reactions */
  while (temp > END) {
      temp = 0;
      for (i = 0; i < nbEspeces; i++) {
          tempo = SBML_simulate(mod, molecules, r, TR, nbEspeces, i);
          /*if (tempo != END) {
              Especes_print_2(molecules, nbEspeces);
				 printf("\n");
          }*/
          temp += tempo;
      }
  }
  printf("Fin de simulation...\n");
  printf("\nEtat des especes :\n");
  Especes_print_2(molecules, nbEspeces);
  /* Liberation de la memoire des generateurs aleatoire */
  gsl_rng_free(r);

  /* Liberation de la memoire de la structure Especes */
  Especes_free(molecules, nbEspeces);
  if (TR != NULL)
    free(TR);
}

/**
 * \fn void SBML_score_add(pScore result, pScore result_temp, FILE *debugFile)
 * \author Amine Ghozlane
 * \brief  Add scores
 * \param  result Struct Score used for all the simulation
 * \param  result_temp Struct Score used at each simulation step
 * \param  debugFile File use for debug
 */
void SBML_score_add(pScore result, pScore result_temp, int nbEspeces)
{
  /* Addition des scores */
  int i;

  /* Copie des resultats */
  for(i=0;i<nbEspeces;i++){
      if(result->name[i]==NULL){
          result->name[i]=(char*)malloc(((int)strlen(result_temp->name[i])+1)*sizeof(char));
          assert(result->name[i]!=NULL);
          strcpy(result->name[i], result_temp->name[i]);
      }
      result->quantite[i]+=result_temp->quantite[i];
  }
}

/**
 * \fn void SBML_score_mean(pScore result, int n)
 * \author Amine Ghozlane
 * \brief  Mean quantities for score
 * \param  result Struct Score
 * \param  n Number of simulation step
 */
void SBML_score_mean(pScore result, int nbEspeces, int n)
{
  /* Moyenne des resultats */
  int i;
  for(i=0;i<nbEspeces;i++){
      result->quantite[i]/=n;
  }
}

/**
 * \fn pScore SBML_scoreAlloc(pListParameters a)
 * \author Amine Ghozlane
 * \brief  Allocation of the struct Score
 * \param  a Global parameters : struct ListParameters
 * \return Allocated struct Score
 */
pScore SBML_scoreAlloc(int nbEspeces)
{
  int i;
  /* Allocation de la structure de score */
  pScore result=NULL;
  result=(pScore)malloc(1*sizeof(Score));
  assert(result!=NULL);

  result->name=NULL;
  result->quantite=NULL;
  result->name=(char**)malloc(nbEspeces*sizeof(char*));
  assert(result->name!=NULL);
  result->quantite=(double*)malloc(nbEspeces*sizeof(double));
  assert(result->quantite!=NULL);
  for(i=0;i<nbEspeces;i++){
      result->name[i]=NULL;
      result->quantite[i]=0.0;
  }
  return result;
}

/**
 * \fn void SBML_scoreFree(pScore out)
 * \author Amine Ghozlane
 * \brief  Free the struct Score
 * \param  out Struct score
 */
void SBML_scoreFree(pScore out, int nbEspeces)
{
  int i;
  /* Libere la memoire alloue a la structure de score */
  for(i=0;i<nbEspeces;i++){
      if(out->name[i]!=NULL) free(out->name[i]);
  }
  if(out->name!=NULL) free(out->name);
  if(out->quantite!=NULL) free(out->quantite);
  if(out!=NULL) free(out);
}

/* Affiche le score moyen */

/**
 * \fn void SBML_scoreFree(pScore out)
 * \author Amine Ghozlane
 * \brief  Free the struct Score
 * \param  out Struct score
 * \param  nbEspeces
 */
void SBML_score_print(pScore result, int nbEspeces)
{
  int i;
  /* Print head*/
  /*for(i=0;i<nbEspeces;i++){
      printf("%s;",result->name[i]);
  }
  printf("\n");*/
  /* Print Value */
  for(i=0;i<nbEspeces;i++){
      /*printf("Molecule: %s - Biomass: %.0f\n",result->name[i],result->quantite[i]);*/
      printf("%.0f;",result->quantite[i]);
  }
  printf("\n");
}

/**
 * \fn void SBML_compute_simulation_mean(Model_t *mod, int nb_simulation)
 * \author Amine Ghozlane
 * \brief  X time simulation of metabolic network
 * \param  mod Model of the SBML file
 * \param  nb_simulation Number of simulation step
 */
void SBML_compute_simulation_mean(Model_t *mod, int nb_simulation)
{
  /* Simulation du reseau metabolique  */
  int i, j,nbReactions = 0, nbEspeces = 0, temp = 1, tempo=0, test=0;
  pEspeces molecules=NULL;
  pTestReaction TR=NULL;
  const gsl_rng_type *T;
  gsl_rng * r;
  pScore result=NULL,result_temp=NULL;

  /* Allocation memoire et initialisation des generateurs de nombre aleatoire */
  gsl_rng_env_setup();
  if (!getenv("GSL_RNG_SEED"))
    gsl_rng_default_seed = time(NULL)+(double)getpid();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*gsl_rng_default_seed = time(NULL)+(double)getpid();*/

  /* Allocation memoire */
  TR=(pTestReaction)malloc(1*sizeof(TestReaction));
  assert(TR!=NULL);

  /* File information */
  nbReactions = (int)Model_getNumReactions(mod);
  nbEspeces = (int)Model_getNumSpecies(mod);
  molecules = Especes_alloc(nbEspeces);

  /* Initialisation de la structure de score temporaire */
  result=SBML_scoreAlloc(nbEspeces);
  result_temp=SBML_scoreAlloc(nbEspeces);

  /* Initialisation de la quantite des especes */
  SBML_initEspeceAmounts(mod, molecules, nbEspeces);
  /* Initialisation des reactions et des ratios*/
  SBML_setReactions(mod, molecules, nbReactions, nbEspeces);

  /*printf("Start the simulation\n");*/
  /* SIMULATION */
  for(j=0;j<nb_simulation;j++){
      /* Simulation des reactions */
      while (temp > END) {
          temp = 0;
          for (i = 0; i < nbEspeces; i++) {
              tempo = SBML_simulate(mod, molecules, r, TR, nbEspeces, i);
              temp +=tempo;
          }
          test+=1;
      }
      temp=1;
      tempo=0;
      /*printf("Nombre de tour %d\n",test);*/
      /*Score */
      /* Enregistre le score des especes */
      Especes_scoreSpecies(molecules, nbEspeces, result_temp->name, result_temp->quantite);
      SBML_score_add(result,result_temp,nbEspeces);
      SBML_initEspeceAmounts(mod, molecules, nbEspeces);
  }
  SBML_score_mean(result,nbEspeces,nb_simulation);

  /*printf("End of the simulation...\n");
  printf("Final state of the molecules :\n");*/
  SBML_score_print(result, nbEspeces);

  /* Liberation de la memoire des generateurs aleatoire */
  gsl_rng_free(r);

  /* Liberation de la memoire du score */
  SBML_scoreFree(result, nbEspeces);
  SBML_scoreFree(result_temp, nbEspeces);

  /* Liberation de la memoire de la structure Especes */
  Especes_free(molecules, nbEspeces);
  if(TR!=NULL) free(TR);
}
