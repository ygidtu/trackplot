#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/10 2:04 PM
# __author__ = 'Zhou Ran'
u"""
Fetch protein information from uniprot website
"""
from collections import namedtuple
from types import SimpleNamespace

import numpy as np
import requests as rq
import xmltodict
import json
from xml.parsers.expat import ExpatError
from conf.logger import logger


class Uniprot(object):
    u"""
    Get domain information from uniprot based uniprot id
    """

    def __init__(self, uniprot_id: str, cds_len: int, fmt="xml", database="uniprot"):
        u"""
        Fetch transcript's translation information based gene id
        :param uniprot_id: the ensemble transcript id works well currently
        :param cds_len: the length of CDS for checking results from uniprot
        :param fmt: default: xml
        :param database: default: uniprot
        """
        self.ui = uniprot_id
        self.cds_len = cds_len
        self.fmt = fmt
        self.database = database

        self.__valid_fmt = {'txt', 'xml', 'rdf', 'gff', 'fasta'}
        self.__url = "https://www.uniprot.org"
        self.domain = self.__domain_info__()
        self.guessed_id = self.__guess_protein_id__()

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

    def __guess_protein_id__(self):
        u"""
        Use xml to parser the response contents and return a dict
        :return: a dict which contained the response results from the given url.
        """

        try:
            # No Domain information
            if self.__request_url__().text == "":
                return {}

            xml_dic = xmltodict.parse(
                self.__request_url__().text, attr_prefix="", cdata_key=""
            )["uniprot"]["entry"]
            current_uniprot_id = set()

            u"""
            There were more than one protein for a transcript, like ENST00000421241.
            We selected the first protein for downstream visualization because usually was reviewed by uniprot, 
            """

            all_alternative_uniprot_id = []
            if isinstance(xml_dic, list):
                features_length = []
                for sub_domain in xml_dic:
                    features_length.append(
                        int(sub_domain['sequence']['length'])
                    )
                    try:
                        for i in sub_domain["comment"]:
                            if "isoform" in i.keys():
                                all_alternative_uniprot_id.extend(list(map(lambda x: x['id'], i["isoform"])))
                    except KeyError:
                        continue

                features_length = np.array(features_length)
                length_match_index = np.where(
                    features_length == int(self.cds_len / 3)
                )[0]
                u"""
                Some protein's length was not equal to uniprot, the transcript maybe a non-canonical isoform.  
                """
                if len(length_match_index) != 0:
                    # if multiple, then return the value at minimal index
                    all_alternative_uniprot_id = []
                    xml_dic = xml_dic[min(length_match_index)]
                else:
                    xml_dic = xml_dic[0]

            if int(xml_dic["sequence"]["length"]) != int(self.cds_len / 3):

                if len(all_alternative_uniprot_id) == 0:
                    try:
                        for i in xml_dic["comment"]:
                            if "isoform" in i.keys():
                                all_alternative_uniprot_id = list(map(lambda x: x['id'], i["isoform"]))
                    except KeyError:
                        pass

                if len(all_alternative_uniprot_id) == 0:
                    return None

                for alternative_id in all_alternative_uniprot_id:
                    request_url = f"https://www.ebi.ac.uk/proteins/api/features/{alternative_id}"
                    try:
                        current_uniprot_inf = rq.get(
                            request_url,
                            timeout=10
                        )
                        current_uniprot_inf = json.loads(current_uniprot_inf.text)
                        if len(current_uniprot_inf["sequence"]) == self.cds_len / 3:
                            current_uniprot_id.add(current_uniprot_inf['accession'])

                    except ConnectionError:
                        raise f"Timeout for {request_url}."

            else:
                if isinstance(xml_dic['accession'], list):
                    return xml_dic['accession'][0]
                return xml_dic['accession']

        except ExpatError:
            logger.warning(f"Timeout or no domain information found, id: {self.ui}.")
            return None

        # if there were multiple protein id, return first.
        return list(current_uniprot_id)[0]

    def __domain_info__(self):
        u"""
        Check the attribution of "feature" in the response results.
        :return: a list which contained feature's attribution
        """
        feature_info = \
            rq.get(
                f"https://www.ebi.ac.uk/proteins/api/features/{self.__guess_protein_id__()}",
                timeout=10
            )

        try:
            feature_info = json.loads(feature_info.text)['features']
            res = []
            if not isinstance(feature_info, list):
                feature_info = [SimpleNamespace(**feature_info)]

            for sub_feature in feature_info:
                sub_feature = SimpleNamespace(**sub_feature)
                if sub_feature.type != "VAR_SEQ":
                    res.append(sub_feature)

            if len(res) == 0:
                return None

            return res
        except KeyError:
            return None


def main():
    trans_id = 'ENST00000477196'
    trans_id_pep = Uniprot(uniprot_id=trans_id, cds_len=294)

    # trans_id = 'ENST00000379319'
    # trans_id_pep = Uniprot(uniprot_id=trans_id, cds_len=594)

    # trans_id = 'ENST00000486161'
    # # 5736 real, 4515 false
    # trans_id_pep = Uniprot(uniprot_id=trans_id, cds_len=4515)

    print(trans_id_pep.guessed_id)
    print(trans_id_pep.domain)
    for i in trans_id_pep.domain:
        print(i.category, i.type, i.begin, i.end)
    # print(type(trans_id_pep.info))
    # print(trans_id_pep.info)
    # print('done')
    # print(trans_id_pep.ensembl_info, trans_id_pep.__domain_info__)


if __name__ == '__main__':
    main()
