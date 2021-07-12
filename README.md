# LODE
Nella cartella sono presenti alcuni script usati durante la tesi, per lo sviluppo del metodo LODE.

I file LODE_rw e LODE_llt_stg contengono le funzioni per analizzare, col metodo LODE, un qualunque serie storica, modellandola con un random walk o un local linear trend più stagionalità e ricercando valori anomali additivi e cambi di livello. In fondo ai medesimi file è riportato un veloce esempio per il loro utilizzo.

Il file LODE_stg contiene la funzione per analizzare, col metodo LODE, un qualunque serie storica, modellandola con un random walk più stagionalità e ricercando cambi di stagionalità. In fondo al file è riportato un veloce esempio per il suo utilizzo.

I file LODE_Nilo e LODE_UKdeaths analizzano le serie serie storiche della portata annuale d'acqua del Nilo e del numero di guidatori morti o gravemente feriti in UK.

I file sim_e_analisi_rw e sim_e_analisi_llt_stg permettono la riproduzione dei risultati delle simulazioni presentate nella tesi. Rispettivamente per l'analisi di RW con outlier provenienti da misutra di gaussiane e di serie modellate con un local linear trend più stagionalità.

I file sim_2ao_2ls e sim_llt_stg contengono il codice per simulare le serie storiche come descritto nella tesi, nelle sezioni "2 AO e 2 LS" e "Simulazione delle serie storiche".

Il file auxres_Nilo permette di applicare il metodo basato sui residui ausiliati, che abbiamo denominato auxres, descritto nella tesi, alla serie storica del Nilo.
