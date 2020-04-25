# Descrizione
Implementazione struttura dati 2-hop cover tramite tecnica RXL utilizzando l'algoritmo di ordiamento SamPG 

# Librerie necessarie
 - libboost;
 - networkit;
 - python matplotlib.
 
 # Installazione
 
 Modificare il makefile inserendo in INCLUDEPATH i path delle librerie necessarie sul proprio dispositivo e digitare ``` make ``` da terminale.
 
 # Opzioni
 
 - ``` -g ``` path del grafo;
 - ``` -k ``` numero di alberi iniziali;
 - ``` -c ``` numero di contatori;
 - ``` -n ``` numero di nuovi alberi generati ad ogni iterazione;
 - ``` -m ``` numero massimo di alberi;
 - ``` -o ``` path di output;
 - ``` -e ``` modalità di esecuzione {0: run RXL; 
                                      1: test di valori multipli per ciascun iperparam.; 
                                      2: comparazione RXL-PLL; 
                                      3: stampa dei plot;
                                      4: conversione grafi.}
                                      

