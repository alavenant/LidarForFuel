# Installation de l'environnement Conda

Ce fichier décrit comment créer et utiliser l'environnement conda pour le package LidarForFuel.

## Création de l'environnement

Pour créer l'environnement conda à partir du fichier `environment.yml` :

```bash
conda env create -f environment.yml
```

Le fichier `environment.yml` se trouve à la racine du projet.

## Activation de l'environnement

Une fois l'environnement créé, activez-le avec :

```bash
conda activate lidarforfuel
```

## Mise à jour de l'environnement

Si le fichier `environment.yml` a été modifié, vous pouvez mettre à jour l'environnement avec :

```bash
conda env update -f environment.yml --prune
```

## Désactivation de l'environnement

Pour désactiver l'environnement :

```bash
conda deactivate
```

## Suppression de l'environnement

Pour supprimer complètement l'environnement :

```bash
conda env remove -n lidarforfuel
```

## Utilisation

Une fois l'environnement activé, vous pouvez utiliser les modules Python :

```python
# Depuis le dossier python/
import sys
sys.path.append('python')

from fPCpretreatment import fPCpretreatment
from fCBDprofile_fuelmetrics import fCBDprofile_fuelmetrics
from ffuelmetrics import ffuelmetrics
```

Ou si vous êtes dans le dossier `python/` :

```python
from fPCpretreatment import fPCpretreatment
from fCBDprofile_fuelmetrics import fCBDprofile_fuelmetrics
from ffuelmetrics import ffuelmetrics
```

## Notes

- L'environnement utilise `conda-forge` comme canal principal pour la plupart des packages
- Le package `laspy` est installé via pip car il n'est pas toujours disponible dans conda
- `rasterio` est installé via conda depuis conda-forge
- Python 3.8 ou supérieur est requis

