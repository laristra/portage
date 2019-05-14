/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
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

    $('#main-menu li > a[href="index.html"]').text("portage");
    // add icons
    $('#main-menu li > a[href="index.html"]')
	.prepend("<i class='fa fa-cog fa-lg'></i> ");
    $('#main-menu li > a[href="concepts.html"]')
	.prepend("<i class='fa fa-rocket fa-lg'></i> ");
    $('#main-menu li > a[href="example.html"]')
	.prepend("<i class='fa fa-desktop fa-lg'></i> ");
    // these two still use <span> for the dropdown, so leave this.
    $('#main-menu li > a[href="namespaces.html"] > span')
	.before("<i class='fa fa-bars'></i> ");
    // classes
    $('#main-menu li > a[href="annotated.html"] > span')
	.before("<i class='fa fa-book'></i> ");

    // use bootstrap navigation
    $("ul#main-menu").addClass("nav nav-pills nav-fill");
    $("li.current>a").addClass("active");
    $(".sm li").addClass("nav-item");
    $(".sm li>a").addClass("nav-link");

    // Remove a doxygen-created link for keywords in the title
    $("div.header div.title a").contents().unwrap();

    // Tweak the navigation for class/namespace trees
    $("#nav-path").removeClass("navpath");
    $("#nav-path > ul").addClass("breadcrumb");
    $("#nav-path li").addClass("breadcrumb-item");

    // $("table.params").addClass("table");
    // $("div.ingroups").wrapInner("<small></small>");
    $("div.levels").css("margin", "0.5em");
    $("div.levels > span").addClass("btn btn-default btn-xs");
    $("div.levels > span").css("margin-right", "0.25em");

    $("table.directory").addClass("table table-striped");
    $("div.summary > a").addClass("btn btn-default btn-xs");
    // $("table.fieldtable").addClass("table");

    // $(".fragment").addClass("well");
    // $(".memitem").addClass("panel panel-default");
    // $(".memproto").addClass("panel-heading");
    // $(".memdoc").addClass("panel-body");
    // $("span.mlabel").addClass("label label-info");

    // $("table.memberdecls").addClass("table");
    // $("[class^=memitem]").addClass("active");

    // $("div.ah").addClass("btn btn-default");
    // $("span.mlabels").addClass("pull-right");
    // $("table.mlabels").css("width", "100%")
    // $("td.mlabels-right").addClass("pull-right");

    // $("div.ttc").addClass("panel panel-primary");
    // $("div.ttname").addClass("panel-heading");
    // $("div.ttname a").css("color", 'white');
    // $("div.ttdef,div.ttdoc,div.ttdeci").addClass("panel-body");

//    $('#MSearchBox').parent().remove();

    // $('div.fragment.well div.line:first').css('margin-top', '15px');
    // $('div.fragment.well div.line:last').css('margin-bottom', '15px');

    // $('table.doxtable').removeClass('doxtable').addClass('table table-striped table-bordered').each(function(){
    // 	$(this).prepend('<thead></thead>');
    // 	$(this).find('tbody > tr:first').prependTo($(this).find('thead'));

    // 	$(this).find('td > span.success').parent().addClass('success');
    // 	$(this).find('td > span.warning').parent().addClass('warning');
    // 	$(this).find('td > span.danger').parent().addClass('danger');
    // });



    // if($('div.fragment.well div.ttc').length > 0)
    // {
    //     $('div.fragment.well div.line:first').parent().removeClass('fragment well');
    // }

    // $('table.memberdecls').find('.memItemRight').each(function(){
    //     $(this).contents().appendTo($(this).siblings('.memItemLeft'));
    //     $(this).siblings('.memItemLeft').attr('align', 'left');
    // });

    // /*$('table.memberdecls').find('.memTemplItemRight').each(function(){
    //     $(this).contents().appendTo($(this).siblings('.memTemplItemLeft'));
    //     $(this).siblings('.memTemplItemLeft').attr('align', 'left');
    // });*/

    // function getOriginalWidthOfImg(img_element) {
    // 	var t = new Image();
    // 	t.src = (img_element.getAttribute ? img_element.getAttribute("src") : false) || img_element.src;
    // 	return t.width;
    // }

    // $('div.dyncontent').find('img').each(function(){
    // 	if(getOriginalWidthOfImg($(this)[0]) > $('#content>div.container').width())
    // 	    $(this).css('width', '100%');
    // });

    // $(".memitem").removeClass("memitem");
    // $(".memproto").removeClass("memproto");
    // $(".memdoc").removeClass("memdoc");
    // $("span.mlabel").removeClass("mlabel");
    // $("table.memberdecls").removeClass("memberdecls");
    // $("[class^=memitem]").removeClass("memitem");
    // $("span.mlabels").removeClass("mlabels");
    // $("table.mlabels").removeClass("mlabels");
    // $("td.mlabels-right").removeClass("mlabels-right");
    // $(".navpath").removeClass("navpath");
    // $("li.navelem").removeClass("navelem");
    // $("a.el").removeClass("el");
    // $("div.ah").removeClass("ah");
    // $("div.header").removeClass("header");

    // $('.mdescLeft').each(function() {
    // 	if($(this).html()=="&nbsp;") {
    // 	    $(this).siblings('.mdescRight').attr('colspan', 2);
    // 	    $(this).remove();
    // 	}
    // });
    // $('td.memItemLeft').each(function() {
    // 	if($(this).siblings('.memItemRight').html()=="") {
    // 	    $(this).attr('colspan', 2);
    // 	    $(this).siblings('.memItemRight').remove();
    // 	}
    // });

});
    
