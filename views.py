from django.shortcuts import render
from django.http import HttpResponse
from django.template import Context, loader, RequestContext, Template
from django.contrib.auth.decorators import login_required

def index(request):
    context = {}
    if request.user.is_authenticated:
        context["user_authenticated"]=True
        context["username"]=request.user.username
    return render(request, "codon_by_codon/index.html", context)

# This function activates the cgi script.
def calculate(request):
    if request.method == 'POST':
        # Process data a bit
        data = request.POST

        fasta_data = data['functionProtein']
        min_count = data['minCount']
        desc = data['analysisID']
        email_address = data['emailAddress']

        # Run actual calulation (by passing data)
        from . import codon_by_codon
        output_t = codon_by_codon.run(fasta_data, min_count, desc, email_address)
        template = Template(output_t)

        context = RequestContext(request)
        return HttpResponse(template.render(context))
    else:
        return HttpResponse("Please use the form to submit data.")

def help(request):
    context = {}
    if request.user.is_authenticated:
        context["user_authenticated"]=True
        context["username"]=request.user.username
    return render(request, "codon_by_codon/help.html", context)

