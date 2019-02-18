# -*- coding: UTF-8 -*-

"""
Python Code for downloading and checking md5sum
from https://gist.github.com/mjohnsullivan/9322154
"""

import os.path
try:
    import urllib2
except:
    import urllib.request as urllib2
import shutil
import hashlib
import logging
import ssl
from tqdm import tqdm

logger = logging.getLogger('Darts_DNN.get_data')

def validate_file(file_path, hash):
    """
    Validates a file against an MD5 hash value
    :param file_path: path to the file for hash validation
    :type file_path:  string
    :param hash:      expected hash value of the file
    :type hash:       string -- MD5 hash value
    """
    m = hashlib.md5()
    with open(file_path, 'rb') as f:
        while True:
            chunk = f.read(1000 * 1000) # 1MB
            if not chunk:
                break
            m.update(chunk)
    return m.hexdigest() == hash


def download_with_resume(url, file_path, hash=None, timeout=10):
    """
    Performs a HTTP(S) download that can be restarted if prematurely terminated.
    The HTTP server must support byte ranges.
    :param file_path: the path to the file to write to disk
    :type file_path:  string
    :param hash: hash value for file validation
    :type hash:  string (MD5 hash value)
    """
    # don't download if the file exists
    if os.path.exists(file_path):
        return
    block_size = 1000 * 1000 # 1MB
    tmp_file_path = file_path + '.part'
    first_byte = os.path.getsize(tmp_file_path) if os.path.exists(tmp_file_path) else 0
    logger.debug('Starting download at %.1fMB' % (first_byte / 1e6))
    file_size = -1
    gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)  
    try:
        file_size = int(urllib2.urlopen(url).info().get('Content-Length', -1))        
    except:
        file_size = int(urllib2.urlopen(url, context=gcontext).info().get('Content-Length', -1))
    logger.debug('File size is %s' % file_size)
    with tqdm(total=file_size, unit='B', unit_scale=True, unit_divisor=1024) as pbar:
        while first_byte < file_size:
            last_byte = first_byte + block_size \
                if first_byte + block_size < file_size \
                else file_size
            logger.debug('Downloading byte range %d - %d' % (first_byte, last_byte))
            # create the request and set the byte range in the header
            req = urllib2.Request(url)
            req.headers['Range'] = 'bytes=%s-%s' % (first_byte, last_byte)
            try:
                data_chunk = urllib2.urlopen(req, timeout=timeout).read()
            except:
                data_chunk = urllib2.urlopen(req, timeout=timeout, context=gcontext).read()
            # Read the data from the URL and write it to the file
            with open(tmp_file_path, 'ab') as f:
                f.write(data_chunk)
            first_byte = last_byte + 1
            # update pbar
            pbar.update( block_size )

    # rename the temp download file to the correct name if fully downloaded
    if file_size == os.path.getsize(tmp_file_path):
        # if there's a hash value, validate the file
        if hash and not validate_file(tmp_file_path, hash):
            raise Exception('Error validating the file against its MD5 hash')
        shutil.move(tmp_file_path, file_path)
    elif file_size == -1:
        raise Exception('Error getting Content-Length from server: %s' % url)