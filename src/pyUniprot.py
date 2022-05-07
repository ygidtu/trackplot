#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/10 2:04 PM
__author__ = 'Zhou Ran'

from collections import namedtuple
import requests as rq
import xmltodict
from xml.parsers.expat import ExpatError
from conf.logger import logger
from src.GenomicLoci import GenomicLoci


class Uniprot(object):
    u"""
    Get domain information from uniprot based uniprot id
    """

    def __init__(self, uniprot_id: str, fmt="XML", database="uniprot"):
        self.ui = uniprot_id
        self.fmt = fmt
        self.database = database

        self.__valid_fmt = {'txt', 'XML', 'rdf', 'gff', 'fasta'}
        self.__url = "https://www.uniprot.org"
        self.info = self.__info__()

    def __request_url__(self):
        """
        Request the uniprot url and return
        :return: a Response object from requests
        """
        assert self.fmt in self.__valid_fmt, f"Nonlegal format was found, {self.fmt}"

        request_url = f"{self.__url}/uniprot/?query={self.ui}&format={self.fmt}"

        try:
            url_response = rq.get(request_url, timeout=10)
        except ConnectionError:
            raise f"Timeout for {request_url}."
        return url_response

    def __info__(self):
        u"""
        Use xml to parser the response contents and return a dict
        :return: a dict which contained the response results from the given url.
        """
        try:
            xml_dic = xmltodict.parse(
                self.__request_url__().text, attr_prefix="", cdata_key=""
            )["uniprot"]["entry"]
            return xml_dic
        except ExpatError:
            logger.warning(f"Timeout or no domain information found, id: {self.ui}.")
            return None

    @property
    def feature(self):
        u"""
        Check the attribution of "feature" in the response results.
        :return: a list which contained feature's attribution
        """
        try:
            feature_info = self.info['feature']
            res = []
            if not isinstance(feature_info, list):
                res.append(feature_info)
                return res
            return self.info['feature']
        except KeyError:
            return None

    @property
    def __domain_info__(self):
        u"""
        Fetch the domain information
        :return:
        """
        if not self.feature:
            return None
        domain_info = namedtuple('domain_info', ['name', 'type', 'start', 'end'])
        domain_res = []

        for sub_feature in self.feature:
            if 'description' in sub_feature:
                description = 'description'
            else:
                description = 'type'

            current_type = 'type'

            if 'position' in sub_feature['location']:
                domain_res.append(domain_info._make([sub_feature[description],
                                                     sub_feature[current_type],
                                                     int(sub_feature['location']['position']['position']),
                                                     int(sub_feature['location']['position']['position'])]))
            else:
                try:
                    domain_res.append(domain_info._make([sub_feature[description],
                                                         sub_feature[current_type],
                                                         int(sub_feature['location']['begin']['position']),
                                                         int(sub_feature['location']['end']['position'])]))
                except KeyError:
                    '''
                    1.24 there was a location error, because the location status was unknown!!
                    '''
                    if 'status' in sub_feature['location']['begin'] or 'status' in sub_feature['location']['end']:
                        # TODO, here should not return a continue. Return a error class
                        continue

        return domain_res

    @property
    def db_reference(self):
        """
        OrderedDict([('type', 'PROSITE'),
        ('id', 'PS50095'), ('property',
        [OrderedDict([('type', 'entry name'),('value', 'PLAT')]),
        OrderedDict([('type', 'match status'), ('value', '6')])])])
        :return: a nested dict of list
        """
        return self.info['dbReference']

    @property
    def ensembl_info(self):
        """
        OrderedDict([('type', 'Ensembl'), ('id', 'ENSMUST00000208660'),
        ('property', [OrderedDict([('type', 'protein sequence ID'),
        ('value', 'ENSMUSP00000146439')]),
        OrderedDict([('type', 'gene ID'), ('value', 'ENSMUSG00000025900')])])])
        :return: ensemblinfo
        """
        for info in self.db_reference:
            if info['type'] == 'Ensembl':
                return info['id']
        return ""


def main():
    trans_id = 'ENST00000486161'
    trans_id_pep = Uniprot(trans_id)
    print(trans_id_pep.feature)
    print('done')
    print(trans_id_pep.ensembl_info, trans_id_pep.__domain_info__)


if __name__ == '__main__':
    main()
