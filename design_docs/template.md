<h1> Title: Some descriptive title <br>
MPAS-Analysis Team <br>
date: YYYY/MM/DD <br>
</h1>
<h2> Summary </h2>
The purpose of this section is to summarize what capability is to be added to
the MPAS-Analysis system through this design process. It should be clear
what new code will do that the current code does not. Summarizing the primary
challenges with respect to software design and implementation is also
appropriate for this section. Finally, this statement should contain a general
statement with regard to what is "success."


<h1> Requirements </h1>

<h2> Requirement: name-of-requirement-here <br>
Date last modified: YYYY/MM/DD <br>
Contributors: (add your name to this list if it does not appear)
</h2>
Each requirement is to be listed under a "section" heading, as there will be a
one-to-one correspondence between requirements, design, proposed implementation
and testing. Requirements should not discuss technical software issues, but
rather focus on model capability. To the extent possible, requirements should
be relatively independent of each other, thus allowing a clean design solution,
implementation and testing plan.


<h1> Algorithmic Formulations (optional) </h1>
<h2> Design solution: short-description-of-proposed-solution-here <br>
Date last modified: YYYY/MM/DD <br>
Contributors: (add your name to this list if it does not appear)
</h2>
For each requirement, there is a design solution that is intended to meet that
requirement. Design solutions can include detailed technical discussions of
PDEs, algorithms, solvers and similar, as well as technical discussion of
performance issues. In general, this section should steer away from a detailed
discussion of low-level software issues such as variable declarations,
interfaces and sequencing.


<h1> Design and Implementation </h1>
<h2> Implementation: short-desciption-of-implementation-here <br>
Date last modified: YYYY/MM/DD <br>
Contributors: (add your name to this list if it does not appear)
</h2>
This section should detail the plan for implementing the design solution for
requirement XXX. In general, this section is software-centric with a focus on
software implementation. Pseudo code is appropriate in this section. Links to
actual source code are appropriate. Project management items, such as git
branches, timelines and staffing are also appropriate. Pseudo code can be
included via blocks like
```python
def example_function(foo):
    return foo**2.0
```

<h1> Testing </h1>
<h2> Testing and Validation: short-desciption-of-testing-here <br>
Date last modified: YYYY/MM/DD <br>
Contributors: (add your name to this list if it does not appear)
</h2>
How will XXX be tested, i.e., how will be we know when we have met requirement
XXX? What testing will be included for use with `py.test` for continuous integration?
Description of how testing that requires off-line or specialized setup will be used.
