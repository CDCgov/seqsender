from shiny import App, Inputs, Outputs, Session, render, req, ui, reactive
from htmltools import TagList, div

################## INDEX PAGE #######################
index_body = [
    ui.HTML(
        """
<h2>Overview</h2>
<p><code>seqsender</code> is a Python program that is developed to automate the process of generating necessary submission files and batch uploading them to <ins>NCBI archives</ins> (such as <strong>BioSample</strong>, <strong>SRA</strong>, and <strong>Genbank</strong>) and <ins>GISAID databases</ins> (e.g. <strong>EpiFlu</strong>, <strong>EpiPox</strong>, and <strong>EpiCoV</strong>).</p>
<h2>Public Domain Standard Notice</h2>
<p>This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the <a href="https://creativecommons.org/publicdomain/zero/1.0/">CC0 1.0 Universal public domain dedication</a>. All contributions to this repository will be released under the CC0 dedication. By submitting a pull request you are agreeing to comply with this waiver of copyright interest.</p>
<h2>License Standard Notice</h2>
<p>The repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later.</p>
<p>This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version.</p>
<p>This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details.</p>
<p>You should have received a copy of the Apache Software License along with this program. If not, see <a href="http://www.apache.org/licenses/LICENSE-2.0.html">http://www.apache.org/licenses/LICENSE-2.0.html</a></p>
<p>The source code forked from other open source projects will inherit its license.</p>
<h2 id="privacy-standard-notice">Privacy Standard Notice<a href="#privacy-standard-notice"></a>
</h2>
<p>This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the <a href="https://github.com/CDCgov/template/blob/master/DISCLAIMER.md">Disclaimer</a> and <a href="https://github.com/CDCgov/template/blob/master/code-of-conduct.md">Code of Conduct</a>. For more information about CDC’s privacy policy, please visit <a href="https://www.cdc.gov/other/privacy.html">http://www.cdc.gov/other/privacy.html</a>.</p>
<h2 id="contributing-standard-notice">Contributing Standard Notice<a href="#contributing-standard-notice"></a>
</h2>
<p>Anyone is encouraged to contribute to the repository by <a href="https://help.github.com/articles/fork-a-repo">forking</a> and submitting a pull request. (If you are new to GitHub, you might start with a <a href="https://help.github.com/articles/set-up-git">basic tutorial</a>.) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the <a href="http://www.apache.org/licenses/LICENSE-2.0.html">Apache Software License v2</a> or later.</p>
<p>All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at <a href="http://www.cdc.gov/other/privacy.html">http://www.cdc.gov/other/privacy.html</a>.</p>
<h2 id="records-management-standard-notice">Records Management Standard Notice<a href="#records-management-standard-notice"></a>
</h2>
<p>This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the <a href="http://www.cdc.gov">CDC web site</a>.</p>
<h2 id="additional-standard-notices">Additional Standard Notices<a href="#additional-standard-notices"></a>
</h2>
<p>Please refer to <a href="https://github.com/CDCgov/template">CDC’s Template Repository</a> for more information about <a href="https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md">contributing to this repository</a>, <a href="https://github.com/CDCgov/template/blob/master/DISCLAIMER.md">public domain notices and disclaimers</a>, and <a href="https://github.com/CDCgov/template/blob/master/code-of-conduct.md">code of conduct</a>.</p>
"""
    ),
]
