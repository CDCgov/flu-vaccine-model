#### Overview


This repository contains code required to reproduce results from the following paper: 

SE Morris, LA Grohskopf, JM Ferdinands, C Reed, M Biggerstaff (2023) **Evaluating potential impacts of a preferential vaccine recommendation for adults aged 65 and older on United States influenza burden**. Epidemiology.

There is no original data associated with this work. 

*The findings and conclusions of this report are those of the authors and do not necessarily represent the official position of the Centers for Disease Control and Prevention.*

<br>

#### Quick-start guide

The numerical naming of .R files corresponds approximately to the order in which those files are run: 

* Scripts numbered 0--2 construct the basic functions and parameter inputs required to simulate the model (and are used in all subsequent analyses);
* Scripts numbered 3a, 3b simulate the model with baseline parameter values;
* Scripts numbered 4a--4c simulate the main results with different assumptions about vaccine effects (indirect or not) and parameter sampling distributions (uniform or not);
* Scripts numbered 5a--5d conduct sensitivity analyses (one-way, two-way and multi-way (via partial-rank correlation coefficient analysis));
* Scripts numbered 6a, 6b create the figures.


<br>

#### Additional information

Questions can be directed to Sinead Morris (run7@cdc.gov) or Matthew Biggerstaff (zmo2@cdc.gov).

<br>

<small>

#### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal public domain dedication. All contributions to this repository will be released under the CC0 dedication. By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

#### License

The project utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later.

This program is free software: you can redistribute it and/or modify it under the terms of the Apache Software License v2, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html.

##### Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the Disclaimer and Code of Conduct. For more information about CDC's privacy policy, please visit http://www.cdc.gov/other/privacy.html.

##### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by forking and submitting a pull request. (If you are new to GitHub, you might start with a basic tutorial.) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the Apache Software License v2 or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at http://www.cdc.gov/other/privacy.html.

##### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the CDC web site.

</small>
