#!/usr/bin/env python3
#
# This confidential and proprietary software may be used only as
# authorised by a licensing agreement from Arm Limited.
#    Copyright 2020 Arm Ltd. All Rights Reserved.
# The entire notice above must be reproduced on all authorised
# copies and copies may only be made to the extent permitted
# by a licensing agreement from Arm Limited.
#

import argparse
import json
import os
import requests

"""
Sign an executable using the authenticode signing service.
See https://confluence.arm.com/display/PE/Code+Signing+User+Guide for details.
"""

AUTHENTICAL_URL = 'http://authenticode.euhpc.arm.com:8087'

def get_signing_request_url(signing_details):
    """
    Send a request to the code signing server with the metadata for our signing job.
    The server will respond with a one-time url we can use to upload the executable which
    this function returns.
    """
    headers = {
        'Content-Type': 'application/json',
    }
    response = requests.post(AUTHENTICAL_URL + '/signed-items', json=signing_details, headers=headers)
    response.raise_for_status()
    return json.loads(response.text)['upload-url']

def sign_binary(url, input_file, output_file):
    """
    Upload an executable to the code signing server and save the response.
    """
    headers = {
        'Content-Type': 'application/octet-stream',
    }
    with open(input_file, 'rb') as f:
        response = requests.put(url, data=f.read(), headers=headers)
    response.raise_for_status()
    out_dir = os.path.dirname(output_file)
    os.makedirs(out_dir, exist_ok=True)
    with open(output_file, 'wb') as f:
        f.write(response.content)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',
                        help='Path to executable to be signed')
    parser.add_argument('output_file',
                        help='Path where signed executable will be written')
    parser.add_argument('--username',
                        default='mali-tools',
                        help='Username used to connect to the server')
    parser.add_argument('--description',
                        default='Mali Offline Compiler',
                        help='Description of the signed content')
    parser.add_argument('--digest-algorithm',
                        default='sha2',
                        help='Algorithm to use when signing the executable')

    args = parser.parse_args()

    signing_details = {
        'username': args.username,
        'description': args.description,
        'cross-sign': True,
        'digest-algorithm': args.digest_algorithm,
    }
    signing_request_url = get_signing_request_url(signing_details)
    sign_binary(signing_request_url, args.input_file, args.output_file)

