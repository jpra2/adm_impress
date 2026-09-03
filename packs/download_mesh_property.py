from __future__ import annotations

import argparse
import hashlib
import json
from pdb import Pdb
import sys
import tarfile
from pathlib import Path

import gdown

from packs.defpaths import BASE_DIR, CACHE_DIR, MANIFEST_FILE, MESH_PROPERTY_DIR


# ============================================================
# Diretórios do projeto
# ============================================================

PROJECT_ROOT = BASE_DIR

# ============================================================
# Configuração
# ============================================================

HASH_CHUNK_SIZE = 1024 * 1024  # 1 MiB
MESH_DIR = MESH_PROPERTY_DIR


# ============================================================
# Manifest
# ============================================================

def load_manifest() -> dict:
    """
    Carrega o mesh_manifest.json.
    """

    if not MANIFEST_FILE.exists():
        raise FileNotFoundError(
            f"Manifest não encontrado: {MANIFEST_FILE}"
        )

    with MANIFEST_FILE.open(
        "r",
        encoding="utf-8",
    ) as file:
        return json.load(file)


# ============================================================
# SHA-256
# ============================================================

def calculate_sha256(path: Path) -> str:
    """
    Calcula o SHA-256 de um arquivo.
    """

    sha256 = hashlib.sha256()

    with path.open("rb") as file:
        while True:

            chunk = file.read(HASH_CHUNK_SIZE)

            if not chunk:
                break

            sha256.update(chunk)

    return sha256.hexdigest()


# ============================================================
# Verificação SHA-256
# ============================================================

def verify_sha256(
    path: Path,
    expected_sha256: str,
) -> bool:
    """
    Verifica se o SHA-256 do arquivo corresponde
    ao valor esperado.
    """

    if not path.exists():
        return False

    print(
        f"Verificando SHA-256 de:\n"
        f"  {path}"
    )

    actual_sha256 = calculate_sha256(path)

    if actual_sha256.lower() == expected_sha256.lower():

        print("SHA-256: OK")

        return True

    print("SHA-256: INCORRETO")

    print(
        f"Esperado:\n"
        f"  {expected_sha256}"
    )

    print(
        f"Obtido:\n"
        f"  {actual_sha256}"
    )

    return False


# ============================================================
# Verificação de tamanho
# ============================================================

def verify_size(
    path: Path,
    expected_size: int | None,
) -> bool:
    """
    Verifica o tamanho do arquivo, caso o manifest
    contenha essa informação.
    """

    if expected_size is None:
        return True

    if not path.exists():
        return False

    actual_size = path.stat().st_size

    if actual_size == expected_size:
        return True

    print(
        f"Tamanho incorreto para {path.name}."
    )

    print(
        f"Esperado: {expected_size} bytes"
    )

    print(
        f"Obtido:    {actual_size} bytes"
    )

    return False


# ============================================================
# Verificação completa de arquivo
# ============================================================

def verify_file(
    path: Path,
    expected_size: int | None,
    expected_sha256: str,
) -> bool:
    """
    Verifica existência, tamanho e SHA-256.
    """

    if not path.exists():
        return False

    if not verify_size(
        path,
        expected_size,
    ):
        return False

    return verify_sha256(
        path,
        expected_sha256,
    )


# ============================================================
# Caminhos
# ============================================================

def get_mesh_path(
    filename: str,
) -> Path:
    """
    Retorna o caminho da malha dentro de mesh/.

    Também impede caminhos maliciosos no manifest.
    """
    
    MESH_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    mesh_root = MESH_DIR.resolve()

    path = (
        MESH_DIR / filename
    ).resolve()

    if not path.is_relative_to(mesh_root):
        raise ValueError(
            f"Caminho de malha invalido: {filename}"
        )

    return path


def get_archive_path(
    filename: str,
) -> Path:
    """
    Retorna o caminho do .tar.gz dentro de mesh_cache/.
    """

    CACHE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    cache_root = CACHE_DIR.resolve()

    path = (
        CACHE_DIR / filename
    ).resolve()

    if not path.is_relative_to(cache_root):
        raise ValueError(
            f"Caminho de archive inválido: {filename}"
        )

    return path


# ============================================================
# Google Drive
# ============================================================

def download_from_google_drive(
    file_id: str,
    destination: Path,
) -> None:
    """
    Baixa um arquivo público do Google Drive.

    O download é feito primeiro para um arquivo .part.
    Só depois de concluído ele recebe o nome definitivo.
    """

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary_file = destination.with_name(
        destination.name + ".part"
    )

    # Remove eventual download incompleto.
    temporary_file.unlink(
        missing_ok=True
    )

    print()
    print("Baixando do Google Drive...")
    print(
        f"Arquivo: {destination.name}"
    )

    try:

        result = gdown.download(
            id=file_id,
            output=str(temporary_file),
            quiet=False
        )

        if result is None:
            raise RuntimeError(
                "O Google Drive não permitiu "
                "o download."
            )

        temporary_file.replace(
            destination
        )

    except Exception:

        temporary_file.unlink(
            missing_ok=True
        )

        raise


# ============================================================
# Extração segura
# ============================================================


def extract_mesh(
    archive_path: Path,
    expected_filename: str,
) -> Path:
    """
    Extrai a malha esperada do arquivo .tar.gz.

    A malha é extraída para a pasta mesh/.
    Apenas o arquivo esperado é extraído.
    """

    MESH_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    mesh_root = MESH_DIR.resolve()

    print()
    print(
        f"Extraindo {archive_path.name}..."
    )

    with tarfile.open(
        archive_path,
        mode="r:gz",
    ) as archive:

        # Procura a malha dentro do archive
        expected_member = None

        for member in archive.getmembers():

            if member.isdir():
                continue

            if Path(member.name).name == expected_filename:

                if expected_member is not None:
                    raise RuntimeError(
                        f"O arquivo '{expected_filename}' "
                        "aparece mais de uma vez no archive."
                    )

                expected_member = member

        if expected_member is None:
            raise RuntimeError(
                f"A malha '{expected_filename}' "
                f"não foi encontrada em "
                f"'{archive_path.name}'."
            )

        # Caminho de destino
        target = (
            MESH_DIR / expected_filename
        ).resolve()

        # Proteção contra path traversal
        if not target.is_relative_to(mesh_root):
            raise RuntimeError(
                "Caminho de destino inválido."
            )

        # Arquivo temporário
        temporary_file = target.with_name(
            target.name + ".part"
        )

        temporary_file.unlink(
            missing_ok=True
        )

        # Abre o arquivo dentro do tar
        source = archive.extractfile(
            expected_member
        )

        if source is None:
            raise RuntimeError(
                f"Não foi possível extrair "
                f"'{expected_filename}'."
            )

        try:

            # Abre o arquivo temporário
            with temporary_file.open("wb") as output:

                while True:

                    chunk = source.read(
                        HASH_CHUNK_SIZE
                    )

                    if not chunk:
                        break

                    output.write(chunk)

            # Só depois que terminou completamente,
            # transforma o .part no arquivo definitivo.
            temporary_file.replace(
                target
            )

        except Exception:

            temporary_file.unlink(
                missing_ok=True
            )

            raise

        finally:

            source.close()

    print(
        f"Malha extraída para:\n"
        f"  {target}"
    )
    
    archive_path.unlink(missing_ok=True)

    return target


# ============================================================
# Backend
# ============================================================

def download_archive(
    mesh_info: dict,
    destination: Path,
) -> None:
    """
    Faz o download do archive utilizando o backend
    especificado no manifest.
    """

    backend = mesh_info.get(
        "backend"
    )

    if backend == "google_drive":

        file_id = mesh_info.get(
            "file_id"
        )

        if not file_id:
            raise ValueError(
                "file_id não definido no manifest."
            )

        download_from_google_drive(
            file_id=file_id,
            destination=destination,
        )

    elif backend == "zenodo":

        raise NotImplementedError(
            "Backend Zenodo ainda não está "
            "implementado."
        )

    else:

        raise ValueError(
            f"Backend desconhecido: {backend}"
        )


# ============================================================
# Função principal
# ============================================================

def ensure_mesh_property(
    mesh_property_name: str,
) -> Path:
    """
    Garante que a malha esteja disponível localmente.

    Fluxo:

    1. Verifica se .msh existe.
    2. Se existe, verifica SHA-256.
    3. Se estiver correto, retorna.
    4. Se não existe ou está inválido:
       verifica .tar.gz.
    5. Se .tar.gz não existe, baixa.
    6. Verifica SHA-256 do .tar.gz.
    7. Extrai a .msh.
    8. Verifica SHA-256 da .msh.
    9. Retorna o caminho da malha.
    """

    manifest = load_manifest()

    meshes_property = manifest.get(
        "mesh_property",
        {},
    )

    # --------------------------------------------------------
    # Procura a malha no manifest
    # --------------------------------------------------------

    if mesh_property_name not in meshes_property:

        available = "\n".join(
            f"  - {name}"
            for name in meshes_property
        )

        raise ValueError(
            f"Mesh property '{mesh_property_name}' não está "
            f"no manifest.\n\n"
            f"Malhas disponíveis:\n"
            f"{available}"
        )

    mesh_info = meshes_property[mesh_property_name]

    filename = mesh_info["filename"]
    archive_name = mesh_info["archive"]

    mesh_sha256 = mesh_info["mesh_sha256"]
    archive_sha256 = mesh_info["archive_sha256"]

    mesh_size = mesh_info.get(
        "mesh_size"
    )

    archive_size = mesh_info.get(
        "archive_size"
    )

    mesh_path = get_mesh_path(
        filename
    )

    archive_path = get_archive_path(
        archive_name
    )

    # --------------------------------------------------------
    # Informações
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print(
        f"Preparando mesh_propery: {mesh_property_name}"
    )
    print("=" * 70)

    # ========================================================
    # 1. A .msh já existe?
    # ========================================================

    if mesh_path.exists():

        print()
        print(
            f"Mesh property {mesh_property_name} já existe localmente."
        )

        print(
            f"Arquivo:\n"
            f"  {mesh_path}"
        )

        print()
        print(
            "Verificando integridade da malha..."
        )

        mesh_ok = verify_file(
            path=mesh_path,
            expected_size=mesh_size,
            expected_sha256=mesh_sha256,
        )

        if mesh_ok:

            print()
            print(
                "Malha válida."
            )

            print(
                "Nenhum download necessário."
            )

            return mesh_path

        # ----------------------------------------------------
        # .msh existe, mas está inválida
        # ----------------------------------------------------

        print()
        print(
            "A malha local está inválida."
        )

        print(
            "Ela será removida e reconstruída "
            "a partir do archive."
        )

        mesh_path.unlink(
            missing_ok=True
        )

    # ========================================================
    # 2. A .tar.gz existe?
    # ========================================================

    if archive_path.exists():

        print()
        print(
            "O arquivo .tar.gz já existe localmente."
        )

        print(
            f"Arquivo:\n"
            f"  {archive_path}"
        )

        print()
        print(
            "Verificando integridade do archive..."
        )

        archive_ok = verify_file(
            path=archive_path,
            expected_size=archive_size,
            expected_sha256=archive_sha256,
        )

        if not archive_ok:

            print()
            print(
                "O .tar.gz está inválido."
            )

            print(
                "Ele será removido e baixado novamente."
            )

            archive_path.unlink(
                missing_ok=True
            )

    # ========================================================
    # 3. Se não existe, baixar
    # ========================================================

    if not archive_path.exists():

        print()
        print(
            "O .tar.gz não está disponível localmente."
        )

        print(
            "Iniciando download..."
        )

        download_archive(
            mesh_info=mesh_info,
            destination=archive_path,
        )

        # ----------------------------------------------------
        # Verifica o download
        # ----------------------------------------------------

        print()
        print(
            "Verificando o arquivo baixado..."
        )

        archive_ok = verify_file(
            path=archive_path,
            expected_size=archive_size,
            expected_sha256=archive_sha256,
        )

        if not archive_ok:

            archive_path.unlink(
                missing_ok=True
            )

            raise RuntimeError(
                "O .tar.gz baixado não passou "
                "na verificação de integridade."
            )

    # ========================================================
    # 4. O .tar.gz está válido
    # ========================================================

    print()
    print(
        "Archive válido."
    )

    # ========================================================
    # 5. Extrair a malha
    # ========================================================

    mesh_path = extract_mesh(
        archive_path=archive_path,
        expected_filename=filename,
    )

    # ========================================================
    # 6. Verificar a .msh extraída
    # ========================================================

    print()
    print(
        "Verificando a malha extraída..."
    )

    mesh_ok = verify_file(
        path=mesh_path,
        expected_size=mesh_size,
        expected_sha256=mesh_sha256,
    )

    if not mesh_ok:

        mesh_path.unlink(
            missing_ok=True
        )

        raise RuntimeError(
            "A malha extraída não passou "
            "na verificação de integridade."
        )

    # ========================================================
    # 7. Tudo certo
    # ========================================================

    print()
    print("=" * 70)
    print(
        "Mesh Propery pronto"
    )
    print("=" * 70)

    print(
        f"Arquivo:\n"
        f"  {mesh_path}"
    )

    return mesh_path


# ============================================================
# CLI
# ============================================================

def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Gerenciador de malhas do projeto."
        )
    )

    parser.add_argument(
        "mesh",
        help=(
            "Nome da malha definido "
            "no mesh_manifest.json."
        ),
    )

    args = parser.parse_args()

    try:

        mesh_path = ensure_mesh_property(
            args.mesh
        )

        print()
        print(
            f"Mesh property disponivel em:\n"
            f"{mesh_path}"
        )

    except Exception as error:

        print(
            f"\nERRO: {error}",
            file=sys.stderr,
        )

        sys.exit(1)


if __name__ == "__main__":
    main()