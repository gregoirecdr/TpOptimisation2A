# Projet d'optimisation numérique

## Récupérer les sources

Pour récupérer les sources, il faut cloner ce dépot git. Vous pouvez suivre les instructions :
[Cloner un dépôt Gitlab](https://gitlab.irit.fr/toc/etu-n7/documentation/-/wikis/Cloner-un-projet-Gitlab).

## Organisation des sources

Les sujets liés au projet se trouvent dans les notebooks du répertoire `notebooks`. Voici l'ordre des sujets :

* Newton
* Régions de confiance
* Lagrangien augmenté

Pour réaliser le projet vous aurez besoin de cette [documentation](doc-projet.pdf) qui décrit les différents algorithmes à implémenter.

## Rendu

Les consignes complètes et la date limite sont dans la [documentation](doc-projet.pdf). En résumé, lorsque vous avez terminé :

1. exécutez les trois notebooks de bout en bout (*Run All*) puis enregistrez-les : ils doivent être rendus **exécutés**, avec leurs sorties
   (résultats, figures) et vos réponses aux questions ;
2. renommez le répertoire cloné (`projet-optinum`) en `nom_prenom_groupe` ;
3. compressez-le avec `tar -cvzf nom_prenom_groupe.tgz nom_prenom_groupe` et déposez le fichier `.tgz` dans Moodle (les fichiers `.rar` sont à prohiber).

## Installation de Julia

Si vous n'avez pas encore installer Julia sur votre compte N7 ou sur votre machine personnelle, suivez les instructions :
[Utilisation de Julia à l'N7](https://gitlab.irit.fr/toc/etu-n7/documentation/-/wikis/Utilisation-de-Julia-à-l'N7).

### Versions et paquets

Le code et les tests ne nécessitent que les bibliothèques standard `LinearAlgebra` et `Test` ; le paquet `Plots`
est utilisé dans les notebooks pour les graphiques et `IJulia` pour les ouvrir. Ils sont disponibles dans
l'environnement Julia 1.10 de l'N7. Le projet a été testé avec Julia 1.10, 1.12 et 1.13.

Les notebooks s'ouvrent avec le noyau Julia de votre installation : si le noyau proposé a un autre nom (par exemple
`Julia 1.10`), choisissez-le dans le sélecteur de noyau de Jupyter ou de VS Code.

## Illustration de la méthode de Newton

![Newton](newton.gif)
