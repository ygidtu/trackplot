#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import ast
import configparser
import sys

from conf.drawing import RGB
from conf.logger import logger


def default_hive_settings() -> dict:
    """
    Creates dictionary containing the default hive plot settings
    """

    return {'dimension': 1000, 'output_file_path': 'hive.svg', 'draw_hive_plot': True, 'tick_marks': True,
            'axis_subdivision': 0.05, 'tick_labels': True, 'tick_label_font_size': 0, 'tick_label_distance': 3,
            'tick_height': 3, 'tick_thickness': 1, 'axis_start_radius': 50, 'axis_end_radius': 600,
            'axis_thickness': 40, 'axis_angles': [0], 'custom_scale': [], 'draw_bars': True,
            'axis_colors': '#000000', 'bezier_colors': ['#FF0000', '#00FF00', '#0000FF'],
            'use_custom_axis_labels': False, 'axis_labels': None, 'bezier_thickness': 10, 'include_key': True,
            'key_position': ['-700', '-700'], 'key_font_size': 10, 'key_title': 'Key',
            'key_text_color': [0, 0, 0], 'key_title_size': 15}


def default_struct_settings() -> dict:
    """
    Creates dictionary containing the default structure plot settings
    """

    return {'draw_struct_plot': True, 'output_file_path': 'structure.svg', 'plot_width': 1000,
            'plot_height': 600, 'left_margin': 100, 'right_margin': 100, 'top_margin': 100,
            'bottom_margin': 100, 'colors': ['#FF0000', '#00FF00', '#0000FF'], 'axis_color': '#000000',
            'axis_thickness': 2, 'tick_length': 10, 'horiz_label_size': 15, 'horiz_label_spacing': 20,
            'horiz_axis_title_size': 30, 'use_vertical_ticks': True, 'vertical_tick_spacing': 0.2,
            'vert_label_size': 15, 'vert_label_spacing': 20, 'include_key': False, 'key_title_size': 10,
            'key_position': [800, 200], 'key_font_size': 15, 'key_text_color': '#000000'}


def default_sashimi_settings():
    return {
        'width': 7, 'height': 5,
        'intron_scale': 1, 'exon_scale': 1,
        'colors': [
            '#FF8000', '#00C866', '#3399FF', '#db5f57', '#dbc257',
            '#91db57', '#57db80', '#57d3db', '#5770db', '#a157db',
            '#db57b2'
        ],
        'ymax': None, 'number_junctions': True,
        'resolution': 0.5, 'junction_log_base': 10,
        'reverse_minus': False, 'font_size': 6,
        'nyticks': 3, 'nxticks': 4,
        'show_ylabel': True, 'show_xlabel': True,
        'plot_title': None, 'numbering_font_size': 6
    }


def parse_hive_plot_settings(config):
    """
    Parses custom settings for a hive plot from a configuration file. Returns a dictionary containing all the settings

    config is a ConfigParser object which already contains all of the settings
    data is a pandas.DataFrame object containing all of the alternative splicing data
    """
    logger.info('Parsing settings for hive plot...')

    settings = default_hive_settings()

    FLOAT_PARAMS = {
        'dimension',
        'axis_subdivision',
        'tick_label_font_size',
        'tick_label_distance',
        'axis_start_radius',
        'axis_end_radius',
        'bezier_thickness',
        'axis_thickness',
        'tick_height',
        'tick_thickness',
        'axis_label_size',
        'key_font_size',
        'key_title_size'
    }

    BOOLEAN_PARAMS = {
        'tick_marks',
        'tick_labels',
        'draw_bars',
        'include_key',
        'draw_hive_plot'
    }

    OTHER_PARAMS = {
        'axis_colors',
        'bezier_colors',
        'axis_angles',
        'custom_scale',
        'axis_label_radius',
        'key_position',
        'key_text_color',
        'key_position'
    }

    for option in config.options('hive_plot'):
        if option in FLOAT_PARAMS:
            settings[option] = config.getfloat('hive_plot', option)
        elif option in BOOLEAN_PARAMS:
            settings[option] = config.getboolean('hive_plot', option)
        elif option in OTHER_PARAMS:
            settings[option] = ast.literal_eval(config.get('hive_plot', option))

    # check to make sure that the settings have the right format and preprocess them
    colorified = [None] * len(settings['bezier_colors'])
    try:
        for i in range(len(colorified)):
            colorified[i] = RGB.from_hex_string(settings['bezier_colors'][i])
        settings['bezier_colors'] = colorified
    except Exception:
        logger.error('Invalid colors in bezier_colors')
        sys.exit(1)

    try:
        settings['axis_colors'] = RGB.from_hex_string(settings['axis_colors'])
    except Exception:
        logger.error('Invalid color in axis_colors')
        sys.exit(1)

    if settings['custom_scale']:
        if type(settings['custom_scale']) is not list:
            logger.error('custom_scale must be list')
            sys.exit(1)
        for item in settings['custom_scale']:
            if type(item) is not list:
                logger.error('Items in custom_scale must be list')
                sys.exit(1)
            elif len(item) != 2:
                logger.error('Invalid number of elements in element of custom_scale')
                sys.exit(1)
            elif item[0] < 0 or item[0] > 1 or item[1] < item[0] or item[1] > 1:
                logger.error('Invalid boundary in custom_scale')
                sys.exit(1)

    if settings['include_key']:
        if type(settings['key_position']) is not list:
            logger.error('key_position must be list')
            sys.exit(1)
        if len(settings['key_position']) != 2:
            logger.error('key_position can have exactly 2 coordinates')
            sys.exit(1)
        try:
            settings['key_text_color'] = RGB.from_hex_string(settings['key_text_color'])
        except:
            logger.error('key_text_color is not a valid color')

    return settings


def parse_struct_plot_settings(config_parser) -> dict:
    """
    Parses custom settings for a structure plot from a configuration file. Returns a dictionary containing
    all the settings

    config is a ConfigParser object which already contains all of the settings
    data is a pandas.DataFrame object containing all of the alternative splicing data
    """
    logger.info('Parsing settings for structure plot...')

    settings = default_struct_settings()

    FLOAT_PARAMS = {
        'plot_width',
        'plot_height',
        'left_margin',
        'right_margin',
        'top_margin',
        'bottom_margin',
        'axis_thickness',
        'tick_length',
        'horiz_label_size',
        'horiz_label_spacing',
        'horiz_axis_title_size',
        'horiz_axis_title_spacing',
        'vertical_tick_spacing',
        'vert_label_size',
        'vert_label_spacing',
        'key_title_size',
        'key_font_size'
    }

    BOOLEAN_PARAMS = {
        'draw_struct_plot',
        'use_vertical_ticks',
        'include_key',
        'use_custom_key_labels'
    }

    OTHER_PARAMS = {
        'horiz_axis_title',
        'colors',
        'axis_color',
        'output_file_path',
        'key_text_color',
        'key_position'
    }

    for option in config_parser.options('struct_plot'):
        if option in FLOAT_PARAMS:
            settings[option] = config_parser.getfloat('struct_plot', option)
        elif option in BOOLEAN_PARAMS:
            settings[option] = config_parser.getboolean('struct_plot', option)
        elif option in OTHER_PARAMS:
            settings[option] = ast.literal_eval(config_parser.get('struct_plot', option))

    # check to make sure entries in settings have the correct format

    if type(settings['colors']) is not list:
        logger.error('colors must be a list of 3 element lists')
        sys.exit(1)
    else:
        colorified = [None] * len(settings['colors'])
        try:
            for i in range(len(colorified)):
                colorified[i] = RGB.from_hex_string(settings['colors'][i])
            settings['colors'] = colorified
        except Exception:
            logger.error('Invalid colors in colors')
            sys.exit(1)

    try:
        settings['axis_color'] = RGB.from_hex_string(settings['axis_color'])
    except Exception:
        logger.error('Invalid color in axis_color')
        sys.exit(1)

    if settings['include_key']:
        try:
            settings['key_text_color'] = RGB.from_hex_string(settings['key_text_color'])
        except Exception:
            logger.error('Invalid color in key_text_color')
            sys.exit(1)

        if len(settings['key_position']) != 2:
            logger.error('Invalid number of components in key_position')
            sys.exit(1)

        for i in range(2):
            try:
                settings['key_position'][i] = float(settings['key_position'][i])
            except Exception:
                logger.error('Elements in key_position must be numbers')
                sys.exit(1)

    return settings


def parse_sashimi_settings(config_parser) -> dict:
    settings = default_sashimi_settings()

    FLOAT_PARAMS = {
        'width',
        'height',
        'intron_scale',
        'exon_scale',
        # @2018.12.20 discard this setting 'ymax',
        'resolution',
        'junction_log_base',
        'font_size',
        'numbering_font_size'
    }

    INT_PARAMS = {
        'nyticks',
        'nxticks'
    }

    BOOLEAN_PARAMS = {
        'draw_sashimi_plot',
        'number_junctions',
        'reverse_minus',
        'show_ylabel',
        'show_xlabel'
    }

    OTHER_PARAMS = {'colors'}

    for option in config_parser.options('sashimi_plot'):
        if option in INT_PARAMS:
            settings[option] = config_parser.getint('sashimi_plot', option)
        elif option in FLOAT_PARAMS:
            settings[option] = config_parser.getfloat('sashimi_plot', option)
        elif option in BOOLEAN_PARAMS:
            settings[option] = config_parser.getboolean('sashimi_plot', option)
        elif option in OTHER_PARAMS:
            settings[option] = ast.literal_eval(config_parser.get('sashimi_plot', option))

    return settings


def parse_settings(settings_file) -> dict:
    """ Creates multiple dictionaries containing the settings parsed from a settings file.
    Each type of plot has its own settings dictionary.

    settings_file is the name of the text file containing the settings

    Return values:
    data is a pandas.DataFrame object which contains the alternative splicing data
    hive_plot_settings is a dictionary containing the settings for the hive plot
    struct_plot_settings is a dictionary containing the settings for the structure plot
    """
    try:
        config = configparser.ConfigParser()
        logger.info('Reading settings from {0}...'.format(settings_file))
        config.read(settings_file)

        # hive_plot_settings = parse_hive_plot_settings(config)
        # struct_plot_settings = parse_struct_plot_settings(config)

        return parse_sashimi_settings(config)

        # logger.error('Done reading settings.')
        # return hive_plot_settings, struct_plot_settings, sashimi_plot_settings
    except IOError:
        logger.error('{0} is not a valid file path')
        sys.exit(1)


if __name__ == '__main__':
    pass
