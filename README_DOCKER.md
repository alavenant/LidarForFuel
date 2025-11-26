# Docker pour LidarForFuel

Ce document explique comment utiliser Docker pour exécuter le package Python LidarForFuel.

## Image Docker

L'image Docker contient le code Python installé avec toutes les dépendances nécessaires.

## Construction de l'image

```bash
docker build -t lidarforfuel:latest -f Dockerfile .
```

## Utilisation

### Avec Docker Compose (recommandé)

#### Shell interactif
```bash
docker-compose run --rm lidarforfuel
```

### Avec Docker directement

#### Exécuter une commande
```bash
docker run --rm -v $(pwd)/python:/app/python -v $(pwd)/test:/app/test \
  lidarforfuel:latest python -m pytest test/ -v
```

#### Shell interactif
```bash
docker run --rm -it -v $(pwd)/python:/app/python -v $(pwd)/test:/app/test \
  lidarforfuel:latest
```

#### Exécuter fPCpretreatment
```bash
docker run --rm -v $(pwd)/python:/app/python \
  -v $(pwd)/inst:/app/inst \
  -v $(pwd)/Data:/app/Data \
  lidarforfuel:latest \
  python -m python.fPCpretreatment -i /app/inst/extdata/M30_FontBlanche.laz \
  -o /app/Data/output.laz
```

## Exécuter les tests

```bash
# Avec docker-compose
docker-compose run --rm lidarforfuel pytest test/ -v

# Avec docker directement
docker run --rm -v $(pwd)/python:/app/python -v $(pwd)/test:/app/test \
  lidarforfuel:latest pytest test/ -v
```

## Volumes montés

Les volumes suivants sont montés par défaut dans docker-compose :
- `./python` → `/app/python` : Code source Python
- `./test` → `/app/test` : Tests
- `./inst` → `/app/inst` : Données de test
- `./Data` → `/app/Data` : Données (dev seulement)

## Variables d'environnement

- `PYTHONPATH=/app` : Permet d'importer les modules Python depuis n'importe où

## Notes

- L'image utilise micromamba pour une installation rapide de l'environnement conda
- L'environnement conda `lidarforfuel` est créé à partir de `environment.yml`
- Les fichiers R ne sont pas inclus dans l'image (voir `.dockerignore`)
- L'image démarre avec un shell interactif par défaut pour faciliter le développement

