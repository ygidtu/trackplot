#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/9 10:48 AM

u"""
A wrapper for mygene.info, Re-write by Ran zhou at 2022.4.28

including:
    1. get all transcript/protein id based given gene id, and only support ensemble id currently.


"""

from biothings_client import get_client
from conf.logger import logger


class ClientAttributionError(Exception):
    pass


class AttrDict(dict):
    u"""
    Return an attribution-fetched dict.
    """
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class MyGene:
    def __init__(self,
                 input_query: str,
                 fields: str,
                 annot_func: str,):
        u"""
        A mygene query class
        :param input_query: a request query to mygene.info, only used ensembl id currently
        :param fields: default value was "all"
        :param annot_func: category of annotation, one of {"gene", "variant", "chem", "disease", "taxon"}
        """

        self.category = {"gene", "variant", "chem", "disease", "taxon"}

        self.input_query = input_query
        self.annot_func = annot_func
        self.annot_field = fields

        if self.annot_func not in self.category:
            raise ClientAttributionError(
                "The function \"{}\" can't find in \"{}\"".format(self.annot_func, "; ".join(self.category))
            )

        self.__my_client__ = get_client(self.annot_func)
        self.query = self.__query__()

    def __query__(self):
        query_info = self.__my_client__.query(self.input_query, fields=self.annot_field)["hits"]
        return query_info

    @property
    def unique_query(self):

        u"""
         stop the process, and raise a contribution error for these un-unique hits.
        """

        assert len(self.query) == 1, \
            "The query information was not unique, pls check {}".format(self.input_query)

        return self.query[0]

    @property
    def uniprot_id(self) -> str:
        u"""
        Return uniprot id.
        """

        try:
            return self.unique_query["uniprot"]["Swiss-Prot"]
        except KeyError as err:
            logger.warn(f"No uniprot information found! {err}")
            return ""

    @property
    def ensembl_id(self) -> AttrDict:
        """

        :return: all ensemble information
        """

        try:
            return AttrDict(self.unique_query["ensembl"])
        except KeyError as err:
            logger.warn(f"No ensembl information found! {err}")
            return AttrDict()

    @property
    def location(self) -> dict:
        """
        Return a dict which contained a gene's chromosome, genomic coordinates and strand information. like
        {'chr': '6', 'end': 30213427, 'ensemblgene': 'ENSG00000234127', 'start': 30184455, 'strand': -1}

        :return:

        """
        try:
            return AttrDict(self.unique_query["genomic_pos"])
        except ValueError as err:
            u"""
            1.25 sorted the genomic coordinate, and choose the shortest chromosome as candidate.
            """

            logger.warn(f"Multiple genomic hits were found for {self.input_query}, return the first.")
            first_obj_sorted_by_chromosome = sorted(self.unique_query["genomic_pos"], key=lambda x: x["chr"])[0]
            return AttrDict(first_obj_sorted_by_chromosome)


def main():
    current_test = MyGene("ensembl.gene:ENSMUSG00000051951", "all", "gene")
    print(current_test.location)
    print(current_test.ensembl_id)
    print(current_test.uniprot_id)


if __name__ == '__main__':
    main()
