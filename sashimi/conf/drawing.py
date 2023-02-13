#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# all the drawing functions in this file assume that
# the coordinates of the svg file have been transformed to cartesian coordinates
class RGB(object):

    __slots__ = ["red", "green", "blue"]

    def __init__(self, red=255, green=255, blue=255):
        if red > 255 or red < 0 or green > 255 or green < 0 or blue > 255 or blue < 0:
            raise Exception("Invalid color")
        else:
            self.red = int(red)
            self.green = int(green)
            self.blue = int(blue)

    @classmethod
    def from_hex_string(cls, hex_string):
        if len(hex_string) == 7:
            hex_string = hex_string[1:]

        if len(hex_string) == 6:
            red = int(hex_string[0:2], 16)
            green = int(hex_string[2:4], 16)
            blue = int(hex_string[4:], 16)
            return cls(red, green, blue)
        else:
            raise Exception("not a valid hexadecimal color")

    def __str__(self):
        return 'rgb({0}, {1}, {2})'.format(self.red, self.green, self.blue)


def draw_line(svg_file, x1, y1, x2, y2, thickness=10, color=RGB(0, 0, 0)):
    svg_line = '<line x1="{0}" y1="{1}" x2="{2}" y2="{3}" style="stroke: {4}; stroke-width: {5}"/> \n'.format(x1, y1,
                                                                                                              x2, y2,
                                                                                                              color,
                                                                                                              thickness)
    svg_file.write(svg_line)


def draw_line_polar(svg_file, x_origin, y_origin, length, angle, thickness=10, color=RGB(0, 0, 0)):
    svg_line = '<line x1="{0}" y1="{1}" x2="{2}" y2="{3}" style="stroke: {4}; stroke-width: {5}" transform="rotate({6},{0},{1})"/> \n'.format(
        x_origin, y_origin, x_origin + length, y_origin, color, thickness, angle)
    svg_file.write(svg_line)


def draw_bezier(svg_file, x1, y1, x2, y2, controlX, controlY, color=RGB(0, 0, 0), thickness=1):
    svg_bezier = '<path d = "M{0},{1} Q{4},{5} {2},{3}" fill="none" style="stroke: {6}; stroke-width: {7}"/> \n'.format(
        x1, y1, x2, y2, controlX, controlY, color, thickness)
    svg_file.write(svg_bezier)


def draw_text(svg_file, words, x, y, size, angle=0, color=RGB(0, 0, 0)):
    svg_text = '<text x="{0}" y="{1}" style="fill:{2}" font-size="{3}" font-family="sans-serif" transform="scale(1,-1) rotate({5},{0},{1})" text-anchor="middle">{4}</text>\n'.format(
        x, -y, color, size, words, angle)
    svg_file.write(svg_text)


def draw_text_left(svg_file, words, x, y, size, angle=0, color=RGB(0, 0, 0)):
    svg_text = '<text x="{0}" y="{1}" style="fill:{2}" font-size="{3}" font-family="sans-serif" transform="scale(1,-1) rotate({5},{0},{1})" text-anchor="left">{4}</text>\n'.format(
        x, -y, color, size, words, angle)
    svg_file.write(svg_text)


def draw_multiline_text(svg_file, label, x, y, size, color=RGB(0, 0, 0)):
    words_list = label.split('\n')

    svg_file.write(
        '<text y="{0}" style="fill:{1}" font-size="{2}" font-family="sans-serif" transform="scale(1,-1)" text-anchor="middle">'.format(
            -y, color, size))
    for i in range(len(words_list)):
        if i == 0:
            svg_file.write('<tspan x="{0}">{1}</tspan>'.format(x, words_list[i]))
        else:
            svg_file.write('<tspan x="{0}" dy="{1}">{2}</tspan>'.format(x, size, words_list[i]))
    svg_file.write('</text>')


def draw_rectangle(svg_file, x, y, x_dim, y_dim, fill_color):
    svg_rect = '<rect x="{0}" y="{1}" width="{2}" height="{3}" style="fill:{4}; stroke-width:0.01; stroke:{4}"/>\n'.format(
        x, y, x_dim, y_dim, fill_color);
    svg_file.write(svg_rect)
