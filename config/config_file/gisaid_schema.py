{
    'Submission': {
        'required': True,
        'type': 'dict',
        'schema': {
            'NCBI': {
                'required': False,
                'type': 'dict',
                'schema': {
                    'Username': {
                        'required': False,
                        'type': 'string',
                        'regex': '^[a-zA-Z0-9_.+-]+',
                        'nullable': True
                    },
                    'Password': {
                        'required': False,
                        'type': 'string',
                        'regex': '^[a-zA-Z0-9_.+-]+',
                        'nullable': True
                    },
                    'Table2asn': {
                        'required': False,
                        'type': 'boolean',
                        'nullable': True
                    },
                    'Submission_Position': {
                        'required': False,
                        'type': 'integer',
                        'allowed': [1, 2],
                        'nullable': True
                    },
                    'Description': {
                        'required': False,
                        'type': 'dict',
                        'schema': {
                            'Title': {
                                'required': False,
                                'type': 'string',
                                'nullable': True
                            },
                            'Comment': {
                                'required': False,
                                'type': 'string',
                                'nullable': True
                            },
                            'Organization': {
                                'required': False,
                                'type': 'dict',
                                'schema': {
                                    '@role': {
                                        'required': False,
                                        'type': 'string',
                                        'nullable': True
                                    },
                                    '@type': {
                                        'required': False,
                                        'type': 'string',
                                        'nullable': True
                                    },
                                    '@org_id': {
                                        'required': False,
                                        'type': 'integer',
                                        'nullable': True
                                    },
                                    'Name': {
                                        'required': False,
                                        'type': 'string',
                                        'nullable': True
                                    },
                                    'Address': {
                                        'required': False,
                                        'type': 'dict',
                                        'schema':{
                                            'Affil': {
                                                'required': False,
                                                'type': 'string',
                                                'nullable': True
                                            },
                                            'Div': {
                                                'required': False,
                                                'type': 'string',
                                                'nullable': True
                                            },
                                            'Street': {
                                                'required': False,
                                                'type': 'string',
                                                'nullable': True
                                            },
                                            'City': {
                                                'required': False,
                                                'type': 'string',
                                                'nullable': True
                                            },
                                            'Sub': {
                                                'required': False,
                                                'type': 'string',
                                                'nullable': True
                                            },
                                            'Postal_code': {
                                                'required': False,
                                                'type': 'integer',
                                                'nullable': True
                                            },
                                            'Country': {
                                                'required': False,
                                                'type': 'string',
                                                'nullable': True
                                            },
                                            'Email': {
                                                'required': False,
                                                'type': 'string',
                                                'regex': '^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$',
                                                'nullable': True
                                            },
                                            'Phone': {
                                                'required': False,
                                                'type': 'string',
                                                'nullable': True,
                                                'nullable': True
                                            }
                                        }
                                    },
                                    'Submitter': {
                                        'required': False,
                                        'type': 'dict',
                                        'schema': {
                                            '@email': {
                                                'required': False,
                                                'type': 'string',
                                                'regex': '^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$',
                                                'nullable': True
                                            },
                                            '@alt_email': {
                                                'required': False,
                                                'type': 'string',
                                                'dependencies': ['@email'],
                                                'regex': '^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$',
                                                'nullable': True
                                            },
                                            'Name': {
                                                'required': False,
                                                'type': 'dict',
                                                'schema': {
                                                    'First': {
                                                        'required': False,
                                                        'type': 'string',
                                                        'nullable': True
                                                    },
                                                    'Last': {
                                                        'required': False,
                                                        'type': 'string',
                                                        'nullable': True
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            },
            'GISAID': {
                'required': True,
                'type': 'dict',
                'schema': {
                    'Client-Id': {
                        'required': True,
                        'type': 'string'
                    },
                    'Username': {
                        'required': True,
                        'type': 'string',
                        'regex': '^[a-zA-Z0-9_.+-]+'
                    },
                    'Password': {
                        'required': True,
                        'type': 'string',
                        'regex': '^[a-zA-Z0-9_.+-]+'
                    },
                    'Submission_Position': {
                        'required': False,
                        'type': 'integer',
                        'allowed': [1, 2],
                        'nullable': True
                    }
                }
            }
        }
    }
}
