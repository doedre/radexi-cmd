import os
import ycm_core

flags = [
        '-Wall',
        '-Wextra',
        '-Werror',
        '-fexceptions',
        '-DNDEBUG',
        '-std=c99',
        '-xc',
        '-isystem/usr/include/',
        '-Iinc'
        ]

SOURCE_EXTENSIONS = [ '.c', ]

def Settings(**kwargs):
    return {
            'flags': flags,
            'do_cache': True
            }
