# Descrizione
Implementazione struttura dati ```2-hop cover``` tramite tecnica ```RXL``` utilizzando l'algoritmo di ordinamento ```SamPG``` 

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

## Esempi modalità di esecuzione:

- ``` ./RXL -g "uwGraphs/ff-10000.txt.hist" -k 60 -c 16 -n 10 -m 1500 -o "ff10000.txt" -e 0 ``` →  esegue RXL sul grafo   specificato e produce in output nel file "ff10000.txt" le label associate.

- ```./RXL -g "uwGraphs/socfb-Howard90.mtx" -k 80 100 -c 16 -n 10 20 -m 700 -o "Howard90.txt" -e 1 ``` → esegue RXL con tutte le combinazioni dei valori specificati per gli iperparametri e stampa i risultati nel file "Howard90.txt" (in particolare i nuovi vengono memorizzati in coda a quelli già presenti eventualmente nel file specificato, in questo modo si possono accumulare dati relativi a diversi test e poi realizzare i plot relativi richiamando la modalità di esecuzione "-e 3").

- ```./RXL -g "uwGraphs/socfb-Middlebury45.mtx" -k 20 -c 16 -n 5 -m 500 -o "Middlebury45.txt" -e 2 ``` → esegue sia RXL che PLL sul grafo specificato ed esegue un test di comparazione tra i due (produce in output anche i relativi plot).

- ```./RXL -o "result.txt" -e 3``` → produce i plot relativi ai dati del file "result.txt".

- ```./RXL -g "Graph.hist" -e 4``` → converte il grafo nel file "Graph.hist" nel formato richiesto (solo per grafi non diretti e pesati di http://networkrepository.com).   
                                      

N.B. In generale tutti i plot e file di log vengono salvati nella cartella ```LogFiles ```.
