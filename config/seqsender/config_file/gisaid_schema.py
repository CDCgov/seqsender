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
                        'regex': '\s*[\S]+\s*',
                        'nullable': True
                    },
                    'Password': {
                        'required': False,
                        'type': 'string',
                        'regex': '\s*[\S]+\s*',
                        'nullable': True
                    },
                    'Spuid_Namespace': {
                        'required': False,
                        'type': 'string',
                        'regex': '\s*[\S]+\s*',
                        'nullable': True
                    },
                    'BioSample_Package': {
                        'required': False,
                        'type': 'string',
                        'nullable': True
                    },
                    'GenBank_Auto_Remove_Failed_Samples': {
                        'required': False,
                        'type': 'boolean',
                        'nullable': True
                    },
                    'Publication_Title': {
                        'required': False,
                        'type': 'string',
                        'nullable': True
                    },
                    'Publication_Status': {
                        'required': False,
                        'type': 'string',
                        'regex': '(?i)(\W|^)(unpublished|in-press|published)(\W|$)',
                        'nullable': True
                    },
                    'Submission_Position': {
                        'required': False,
                        'type': 'integer',
                        'allowed': [1, 2],
                        'nullable': True
                    },
                    'Specified_Release_Date': {
                        'required': False,
                        'type': 'string',
                        'regex': '((?i)(\W|^)(\d+\s*(days|weeks|months)|\d{4}-\d{2}-\d{2})(\W|$))|(^\s*$)',
                        'nullable': True
                    },
                    'Link_Sample_Between_NCBI_Databases': {
                        'required': False,
                        'type': 'boolean',
                        'nullable': True
                    },
                    'Description': {
                        'required': False,
                        'type': 'dict',
                        'schema': {
                            'Organization': {
                                'required': False,
                                'type': 'dict',
                                'schema': {
                                    'Role': {
                                        'required': False,
                                        'type': 'string',
                                        'nullable': True
                                    },
                                    'Type': {
                                        'required': False,
                                        'type': 'string',
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
                                            'Postal_Code': {
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
                                            'Email': {
                                                'required': False,
                                                'type': 'string',
                                                'regex': '^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$',
                                                'nullable': True
                                            },
                                            'Alt_Email': {
                                                'required': False,
                                                'type': 'string',
                                                'dependencies': ['Email'],
                                                'regex': '(^\s*$)|(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$)',
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
                        'regex': '\s*[\S]+\s*'
                    },
                    'Password': {
                        'required': True,
                        'type': 'string',
                        'regex': '\s*[\S]+\s*'
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
