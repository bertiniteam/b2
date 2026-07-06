{{ fullname | escape | underline }}

{%- if fullname == "bertini" %}
.. Top-level page: it SHOWS the re-exported flat surface, but the canonical cross-reference
   targets live on the submodule pages (where each object is defined), so mark it :no-index:
   to avoid duplicate-target warnings from the re-exports.
.. automodule:: {{ fullname }}
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
{%- else %}
.. automodule:: {{ fullname }}
   :members:
   :undoc-members:
   :show-inheritance:
{%- endif %}

{% block modules %}
{% if modules %}
.. rubric:: Submodules

.. autosummary::
   :toctree:
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
