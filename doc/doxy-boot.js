/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

$( document ).ready(function() {

    /*
      For some reason Doxygen isn't properly setting the active tab for
      user-defined pages.  This fixes that.
    */
    var thisFile = document.location.pathname.split('/').slice(-1)[0];
    $('li > a[href="' + thisFile + '"]').parent().addClass("current");

    $("div.headertitle").addClass("page-header");
    $("div.title").addClass("h1");

    $('li > a[href="index.html"] > span').before("<i class='fa fa-cog fa-lg'></i> ");
    $('li > a[href="index.html"] > span').text("portage");
    $('li > a[href="quickstart.html"] > span').before("<i class='fa fa-rocket fa-lg'></i> ");
    $('li > a[href="simple_mesh.html"] >span').before("<i class='fa fa-desktop fa-lg'></i> ");

    $('li > a[href="namespaces.html"] > span').before("<i class='fa fa-bars'></i> ");
    $('li > a[href="classes.html"] > span').before("<i class='fa fa-book'></i> ");
    $('li > a[href="annotated.html"] > span').before("<i class='fa fa-book'></i> ");

    $('li > a[href="functions_func.html"] > span').before("<i class='fa fa-list'></i> ");
    $('li > a[href="functions_vars.html"] > span').before("<i class='fa fa-list'></i> ");
    $('li > a[href="functions_enum.html"] > span').before("<i class='fa fa-list'></i> ");
    $('li > a[href="functions_eval.html"] > span').before("<i class='fa fa-list'></i> ");
    $('img[src="ftv2ns.png"]').replaceWith('<span class="label label-danger">N</span> ');
    $('img[src="ftv2cl.png"]').replaceWith('<span class="label label-danger">C</span> ');

    $("ul.tablist").addClass("nav nav-pills nav-justified");
    $("ul.tablist").css("margin-top", "0.5em");
    $("ul.tablist").css("margin-bottom", "0.5em");
    $("li.current").addClass("active");
    $("iframe").attr("scrolling", "yes");

    $("#nav-path > ul").addClass("breadcrumb");

    $("table.params").addClass("table");
    $("div.ingroups").wrapInner("<small></small>");
    $("div.levels").css("margin", "0.5em");
    $("div.levels > span").addClass("btn btn-default btn-xs");
    $("div.levels > span").css("margin-right", "0.25em");

    $("table.directory").addClass("table table-striped");
    $("div.summary > a").addClass("btn btn-default btn-xs");
    $("table.fieldtable").addClass("table");
    $(".fragment").addClass("well");
    $(".memitem").addClass("panel panel-default");
    $(".memproto").addClass("panel-heading");
    $(".memdoc").addClass("panel-body");
    $("span.mlabel").addClass("label label-info");

    $("table.memberdecls").addClass("table");
    $("[class^=memitem]").addClass("active");

    $("div.ah").addClass("btn btn-default");
    $("span.mlabels").addClass("pull-right");
    $("table.mlabels").css("width", "100%")
    $("td.mlabels-right").addClass("pull-right");

    $("div.ttc").addClass("panel panel-primary");
    $("div.ttname").addClass("panel-heading");
    $("div.ttname a").css("color", 'white');
    $("div.ttdef,div.ttdoc,div.ttdeci").addClass("panel-body");

    $('#MSearchBox').parent().remove();

    $('div.fragment.well div.line:first').css('margin-top', '15px');
    $('div.fragment.well div.line:last').css('margin-bottom', '15px');

    $('table.doxtable').removeClass('doxtable').addClass('table table-striped table-bordered').each(function(){
	$(this).prepend('<thead></thead>');
	$(this).find('tbody > tr:first').prependTo($(this).find('thead'));

	$(this).find('td > span.success').parent().addClass('success');
	$(this).find('td > span.warning').parent().addClass('warning');
	$(this).find('td > span.danger').parent().addClass('danger');
    });



    if($('div.fragment.well div.ttc').length > 0)
    {
        $('div.fragment.well div.line:first').parent().removeClass('fragment well');
    }

    $('table.memberdecls').find('.memItemRight').each(function(){
        $(this).contents().appendTo($(this).siblings('.memItemLeft'));
        $(this).siblings('.memItemLeft').attr('align', 'left');
    });

    /*$('table.memberdecls').find('.memTemplItemRight').each(function(){
        $(this).contents().appendTo($(this).siblings('.memTemplItemLeft'));
        $(this).siblings('.memTemplItemLeft').attr('align', 'left');
    });*/

    function getOriginalWidthOfImg(img_element) {
	var t = new Image();
	t.src = (img_element.getAttribute ? img_element.getAttribute("src") : false) || img_element.src;
	return t.width;
    }

    $('div.dyncontent').find('img').each(function(){
	if(getOriginalWidthOfImg($(this)[0]) > $('#content>div.container').width())
	    $(this).css('width', '100%');
    });

    $(".memitem").removeClass("memitem");
    $(".memproto").removeClass("memproto");
    $(".memdoc").removeClass("memdoc");
    $("span.mlabel").removeClass("mlabel");
    $("table.memberdecls").removeClass("memberdecls");
    $("[class^=memitem]").removeClass("memitem");
    $("span.mlabels").removeClass("mlabels");
    $("table.mlabels").removeClass("mlabels");
    $("td.mlabels-right").removeClass("mlabels-right");
    $(".navpath").removeClass("navpath");
    $("li.navelem").removeClass("navelem");
    $("a.el").removeClass("el");
    $("div.ah").removeClass("ah");
    $("div.header").removeClass("header");

    $('.mdescLeft').each(function() {
	if($(this).html()=="&nbsp;") {
	    $(this).siblings('.mdescRight').attr('colspan', 2);
	    $(this).remove();
	}
    });
    $('td.memItemLeft').each(function() {
	if($(this).siblings('.memItemRight').html()=="") {
	    $(this).attr('colspan', 2);
	    $(this).siblings('.memItemRight').remove();
	}
    });

});
    
