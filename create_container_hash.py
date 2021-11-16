#!/usr/bin/env python3

import sys

from xml.etree import ElementTree
from xml.etree.ElementTree import Element

DEFAULT_BASE_IMAGE = "bgruening/busybox-bash:0.1"
def main():
    try:
        base_image = sys.argv[1]
    except IndexError:
        base_image = DEFAULT_BASE_IMAGE
    tool = ElementTree.parse("shm_csr.xml").getroot()
    requirements: Element = tool.find("requirements")
    packages = []
    for req in requirements.findall("requirement"):
        if req.get("type") == "package":
            name = req.text
            version = req.get("version")
            package_string = f"{name}={version}"
            packages.append(package_string)
    with open("container_hash.tsv", mode="wt") as container_hash_file:
        container_hash_file.write("#targets\tbase_image\timage_build\n")
        container_hash_file.write(",".join(packages) + f"\t{base_image}\t0\n")


if __name__ == "__main__":
    main()
